function result = compute_flexibility(cfg, ~, ~, preprocessing)
%COMPUTE_FLEXIBILITY Recreate the frozen Sun-phase campaigns and retuning.
%
% PROTOTYPE:
%   result = rpo2ds.compute_flexibility(cfg, database, targets, preprocessing)
%
% DESCRIPTION:
%   Solve independent Venus/Mars cases from archived NLP setups. Apply the
%   fixed retuning recipe in its recorded order and replace a result only
%   with a converged, lower-delta-v candidate. Failed cases retain diagnostic
%   fields without automatic retries; frozen starts are not final solutions.
%
% INPUT:
%   cfg            Configuration from reproduction_config(); uses
%                  rpo_to_deep_space.flexibility, including retuning recipes.
%   database       Unused positional argument retained by the runner interface.
%   targets        Unused positional argument retained by the runner interface.
%   preprocessing  Schema 3.0 output from load_flexibility_initial_guesses;
%                  includes untuned_setup and retuning_setup for each planet.
%
% OUTPUT:
%   result         venus/mars structs with untuned/tuned campaign arrays and
%                  tuning_log, mirrored under data. Phases are in rad, delta-v
%                  in m/s, and flight time in days. Metrics summarize tuned
%                  cases; completed records the campaign finishing, even if
%                  individual cases failed.

% Validate the archived grid/recipes and select serial or parallel case execution.
settings = cfg.rpo_to_deep_space.flexibility;
overall_timer = tic;
rpo2ds.require_solver(settings.solver_backend);
phase_grid = settings.initial_sun_phase_rad(:);
local_validate_preprocessing(preprocessing, settings, phase_grid);
use_parallel = rpo.parallel_available(settings.use_parallel);
if use_parallel
    fprintf('Parallel Computing Toolbox detected; flexibility cases will run in parallel.\n');
else
    fprintf(['Parallel Computing Toolbox is not available. Flexibility cases will ', ...
        'run serially and take longer; install the toolbox for faster execution.\n']);
end
result = struct();

% Solve each planet independently before applying its recorded retuning recipe.
for target_name = ["venus", "mars"]
    fprintf('Running %s Sun-phase flexibility (%d independent cases).\n', ...
        upper(target_name), numel(phase_grid));
    frozen = preprocessing.(target_name);
    untuned = local_campaign(settings, target_name, phase_grid, frozen, use_parallel);
    tuned = untuned;
    tuning_log = local_tuning_log();

    % Evaluate retuning attempts independently, then apply replacements in recipe order.
    if settings.auto_tune_default
        recipe = cfg.rpo_to_deep_space.flexibility.retuning.(target_name);
        candidates = local_retuning_candidates(settings, target_name, ...
            phase_grid, frozen, recipe, use_parallel);
        for attempt = 1:size(recipe, 1)
            index = recipe(attempt, 1);
            requested_phase = recipe(attempt, 2);
            executed_phase = phase_grid(index);
            run_settings = settings;
            run_settings.nodes_per_day = recipe(attempt, 3);
            run_settings.(['rpo_moon_match_tolerance_' char(target_name) '_km_s']) = recipe(attempt, 4);
            old_time_of_flight = tuned.time_of_flight_days(index);
            candidate = candidates{attempt};
            old_value = tuned.total_delta_v_m_s(index);
            replaced = candidate.converged && (isnan(old_value) || candidate.total_delta_v_m_s < old_value);
            if replaced
                tuned = local_replace_row(tuned, candidate, index);
            end
            tuning_log(end + 1) = local_log_entry(attempt, index, requested_phase, ...
                executed_phase, run_settings, target_name, candidate, old_value, ...
                old_time_of_flight, replaced); %#ok<AGROW>
        end
    end
    result.(target_name) = struct('untuned', untuned, 'tuned', tuned, ...
        'tuning_log', tuning_log);
end

% Store both campaigns and summarize the tuned outcomes.
result.data = struct('venus', result.venus, 'mars', result.mars);
result.configuration = settings;
result.metrics = local_metrics(result);
result.status = struct('completed', true, 'failed_sections', strings(0, 1));
result.provenance.historical_producer = "analyze_rpo2ds_initial_geometry.m (independent Sun-phase campaign)";
result.provenance.elapsed_seconds = toc(overall_timer);
end

function campaign = local_campaign(settings, target_name, phase_grid, frozen, use_parallel)
%LOCAL_CAMPAIGN Solve the untuned frozen phase cases in grid order.
%
% PROTOTYPE:
%   campaign = local_campaign(settings, target_name, phase_grid, frozen, use_parallel)
%
% DESCRIPTION:
%   Use the same case solver for serial and parallel execution. Preserve
%   input order while reporting completed-worker progress; failed cases
%   remain recorded without automatic retry.
%
% INPUT:
%   settings      Frozen flexibility solver/iteration/tolerance settings.
%   target_name   venus or mars; selects logging and matching-tolerance fields.
%   phase_grid    Initial solar phases [rad].
%   frozen        Planet-specific archived NLP setup cells.
%   use_parallel  Logical selecting the independent parallel loop.
%
% OUTPUT:
%   campaign      Struct of phase-indexed delta-v [m/s], TOF [days],
%                 convergence flags, solver diagnostics, and error strings.

campaign = local_empty_campaign(numel(phase_grid));
timer = tic;
rows = cell(numel(phase_grid), 1);
untuned_setup = frozen.untuned_setup;
if use_parallel
    completed_count = 0;
    report_stride = max(1, ceil(numel(phase_grid) / 10));
    progress_queue = parallel.pool.DataQueue;
    afterEach(progress_queue, @report_parallel_progress);
    parfor index = 1:numel(phase_grid)
        rows{index} = local_run_point(settings, index, phase_grid(index), ...
            phase_grid(index), untuned_setup{index});
        send(progress_queue, 1);
    end
else
    for index = 1:numel(phase_grid)
        rows{index} = local_run_point(settings, index, phase_grid(index), ...
            phase_grid(index), untuned_setup{index});
        fprintf('  %s flexibility progress: %d/%d (%.1f s elapsed)\n', ...
            upper(target_name), index, numel(phase_grid), toc(timer));
    end
end
for index = 1:numel(rows)
    campaign = local_replace_row(campaign, rows{index}, index);
end
failed_count = nnz(~campaign.converged);
if failed_count > 0
    fprintf('  %s flexibility: %d case(s) did not converge; continuing without retry.\n', ...
        upper(target_name), failed_count);
end

    function report_parallel_progress(~)
    %REPORT_PARALLEL_PROGRESS Count completed workers and print campaign progress.

        completed_count = completed_count + 1;
        if mod(completed_count, report_stride) == 0 || completed_count == numel(phase_grid)
            fprintf('  %s flexibility progress: %d/%d (%.1f s elapsed)\n', ...
                upper(target_name), completed_count, numel(phase_grid), toc(timer));
        end
    end
end

function candidates = local_retuning_candidates(settings, target_name, ...
        phase_grid, frozen, recipe, use_parallel)
%LOCAL_RETUNING_CANDIDATES Solve the archived retuning attempts independently.
%
% PROTOTYPE:
%   candidates = local_retuning_candidates(settings, target_name, ...
%       phase_grid, frozen, recipe, use_parallel)
%
% DESCRIPTION:
%   Use the same case solver for serial and parallel execution. Preserve
%   input order while reporting completed-worker progress; failed cases
%   remain recorded without automatic retry.
%
% INPUT:
%   settings      Frozen flexibility solver/iteration/tolerance settings.
%   target_name   venus or mars; selects logging and matching-tolerance fields.
%   phase_grid    Initial solar phases [rad].
%   frozen        Planet-specific archived NLP setup cells.
%   recipe        Rows [grid index, requested phase, nodes/day, tolerance].
%   use_parallel  Logical selecting the independent parallel loop.
%
% OUTPUT:
%   candidates    Cell per recipe attempt, containing local_run_point output.
%                 Ordered replacement decisions are made later by the caller.

candidates = cell(size(recipe, 1), 1);
timer = tic;
completed_count = 0;
retuning_setup = frozen.retuning_setup;
grid_indices = recipe(:, 1);
requested_phases = recipe(:, 2);
nodes_per_day = recipe(:, 3);
match_tolerances = recipe(:, 4);
executed_phases = phase_grid(grid_indices);
if use_parallel
    progress_queue = parallel.pool.DataQueue;
    afterEach(progress_queue, @report_parallel_progress);
    parfor attempt = 1:size(recipe, 1)
        index = grid_indices(attempt);
        run_settings = settings;
        run_settings.nodes_per_day = nodes_per_day(attempt);
        run_settings.(['rpo_moon_match_tolerance_' char(target_name) '_km_s']) = ...
            match_tolerances(attempt);
        candidates{attempt} = local_run_point(run_settings, index, ...
            executed_phases(attempt), requested_phases(attempt), ...
            retuning_setup{attempt});
        send(progress_queue, 1);
    end
else
    for attempt = 1:size(recipe, 1)
        index = grid_indices(attempt);
        run_settings = settings;
        run_settings.nodes_per_day = nodes_per_day(attempt);
        run_settings.(['rpo_moon_match_tolerance_' char(target_name) '_km_s']) = ...
            match_tolerances(attempt);
        candidates{attempt} = local_run_point(run_settings, index, ...
            executed_phases(attempt), requested_phases(attempt), ...
            retuning_setup{attempt});
        report_parallel_progress([]);
    end
end
failed_count = nnz(cellfun(@(row) ~row.converged, candidates));
if failed_count > 0
    fprintf('  %s retuning: %d attempt(s) did not converge; continuing without retry.\n', ...
        upper(target_name), failed_count);
end

    function report_parallel_progress(~)
    %REPORT_PARALLEL_PROGRESS Count completed workers and print campaign progress.

        completed_count = completed_count + 1;
        fprintf('  %s retuning progress: %d/%d (%.1f s elapsed)\n', ...
            upper(target_name), completed_count, size(recipe, 1), toc(timer));
    end
end

function row = local_run_point(settings, index, phase_rad, requested_phase, frozen_setup)
%LOCAL_RUN_POINT Solve one frozen setup and record failures without retrying the campaign.

row = struct('index', index, 'phase_rad', phase_rad, 'requested_phase_rad', requested_phase, ...
    'total_delta_v_m_s', NaN, 'time_of_flight_days', NaN, 'converged', false, ...
    'exitflag', NaN, 'objective', NaN, 'constraint_residual_inf', NaN, ...
    'error_identifier', "", 'error_message', "");
try
    solution = rpo2ds.solve_frozen_flexibility_case(frozen_setup, settings);
    row.converged = solution.converged;
    row.exitflag = solution.exitflag;
    row.objective = solution.objective;
    row.constraint_residual_inf = solution.constraint_residual_inf;
    if solution.has_solution
        row.total_delta_v_m_s = solution.total_delta_v_m_s;
        row.time_of_flight_days = solution.time_of_flight_days;
    end
catch exception
    row.error_identifier = string(exception.identifier);
    row.error_message = string(exception.message);
end
end

function local_validate_preprocessing(preprocessing, settings, phase_grid)
%LOCAL_VALIDATE_PREPROCESSING Require the archived schema, phase grid, and retuning recipe to match.

% Fail before the campaign if the frozen optimizer setups do not match it.

valid = isstruct(preprocessing) && ...
    all(isfield(preprocessing, {'configuration', 'venus', 'mars'})) && ...
    isfield(preprocessing, 'schema_version') && ...
    preprocessing.schema_version == "3.0" && ...
    isequal(preprocessing.configuration.phase_grid_rad(:), phase_grid(:)) && ...
    isequal(preprocessing.configuration.venus_retuning, settings.retuning.venus) && ...
    isequal(preprocessing.configuration.mars_retuning, settings.retuning.mars);
for target_name = ["venus", "mars"]
    valid = valid && isfield(preprocessing.(target_name), 'untuned_setup') && ...
        isfield(preprocessing.(target_name), 'retuning_setup') && ...
        numel(preprocessing.(target_name).untuned_setup) == numel(phase_grid) && ...
        numel(preprocessing.(target_name).retuning_setup) == ...
        size(settings.retuning.(target_name), 1);
end
if ~valid
    error('rpo2ds:FlexibilityPreprocessingMismatch', ...
        ['The frozen flexibility preprocessing does not match the configured ', ...
        '36-case campaigns and retuning recipes.']);
end
end

function campaign = local_empty_campaign(count)
%LOCAL_EMPTY_CAMPAIGN Allocate phase-indexed metrics with NaN/false failure placeholders.

campaign = struct('index', (1:count).', 'phase_rad', NaN(count, 1), ...
    'requested_phase_rad', NaN(count, 1), 'total_delta_v_m_s', NaN(count, 1), ...
    'time_of_flight_days', NaN(count, 1), 'converged', false(count, 1), ...
    'exitflag', NaN(count, 1), 'objective', NaN(count, 1), ...
    'constraint_residual_inf', NaN(count, 1), 'error_identifier', strings(count, 1), ...
    'error_message', strings(count, 1));
end

function campaign = local_replace_row(campaign, row, index)
%LOCAL_REPLACE_ROW Store one solved row at its original phase-grid index.

campaign.phase_rad(index) = row.phase_rad;
campaign.requested_phase_rad(index) = row.requested_phase_rad;
campaign.total_delta_v_m_s(index) = row.total_delta_v_m_s;
campaign.time_of_flight_days(index) = row.time_of_flight_days;
campaign.converged(index) = row.converged;
campaign.exitflag(index) = row.exitflag;
campaign.objective(index) = row.objective;
campaign.constraint_residual_inf(index) = row.constraint_residual_inf;
campaign.error_identifier(index) = row.error_identifier;
campaign.error_message(index) = row.error_message;
end

function entry = local_log_entry(attempt, index, requested_phase, executed_phase, settings, target_name, row, old_value, old_time_of_flight, replaced)
%LOCAL_LOG_ENTRY Record requested/executed phase, recipe settings, and replacement outcome.

entry = struct('attempt', attempt, 'grid_index', index, ...
    'requested_phase_rad', requested_phase, 'executed_phase_rad', executed_phase, ...
    'nodes_per_day', settings.nodes_per_day, ...
    'match_tolerance_km_s', local_tolerance(settings, target_name), ...
    'converged', row.converged, 'old_delta_v_m_s', old_value, ...
    'new_delta_v_m_s', row.total_delta_v_m_s, 'old_time_of_flight_days', old_time_of_flight, ...
    'new_time_of_flight_days', row.time_of_flight_days, 'replaced', replaced, ...
    'error_identifier', row.error_identifier, 'error_message', row.error_message);
end

function tolerance = local_tolerance(settings, target_name)
%LOCAL_TOLERANCE Read the target-specific lunar matching tolerance [km/s], or return NaN.

field_name = ['rpo_moon_match_tolerance_' char(target_name) '_km_s'];
if ~isfield(settings, field_name)
    tolerance = NaN;
else
    tolerance = settings.(field_name);
end
end

function log = local_tuning_log()
%LOCAL_TUNING_LOG Create the empty ordered retuning-audit struct array.

log = struct('attempt', {}, 'grid_index', {}, 'requested_phase_rad', {}, ...
    'executed_phase_rad', {}, 'nodes_per_day', {}, 'match_tolerance_km_s', {}, ...
    'converged', {}, 'old_delta_v_m_s', {}, 'new_delta_v_m_s', {}, ...
    'old_time_of_flight_days', {}, 'new_time_of_flight_days', {}, ...
    'replaced', {}, 'error_identifier', {}, 'error_message', {});
end

function metrics = local_metrics(result)
%LOCAL_METRICS Summarize tuned campaign ranges, convergence counts, and failed indices.

for target_name = ["venus", "mars"]
    campaign = result.(target_name).tuned;
    finite = campaign.total_delta_v_m_s(isfinite(campaign.total_delta_v_m_s));
    metrics.(target_name) = struct('converged_count', nnz(campaign.converged), ...
        'failed_indices', find(~campaign.converged), ...
        'delta_v_range_m_s', local_range(finite), ...
        'time_of_flight_range_days', local_range(campaign.time_of_flight_days));
end
end

function values = local_range(array)
%LOCAL_RANGE Return finite extrema, or [NaN, NaN] when no finite values remain.

array = array(isfinite(array));
if isempty(array)
    values = [NaN, NaN];
else
    values = [min(array), max(array)];
end
end
