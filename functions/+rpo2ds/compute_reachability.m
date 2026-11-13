function result = compute_reachability(cfg, database, preprocessing)
%COMPUTE_REACHABILITY Reproduce the frozen fixed-budget directional sweep.
%
% PROTOTYPE:
%   result = rpo2ds.compute_reachability(cfg, database, preprocessing)
%
% DESCRIPTION:
%   Temporarily select the historical fmincon backend and solve each direction
%   from its archived patched initial guess. Restore the caller path on exit.
%   The polar sweep is computed; the plotted surface is reconstructed using
%   the historical symmetry assumption, not additional optimized directions.
%   See docs/REPRODUCIBILITY.md (Run 02) for fresh-run agreement limits.
%
% INPUT:
%   cfg            Configuration from reproduction_config(); uses
%                  rpo_to_deep_space.shared and reachability.
%   database       Dense historical or generated phase-aware RPO-Moon table.
%   preprocessing  Output from load_reachability_initial_guesses; contains
%                  rpo_state_2d and the ordered points matching the polar grid.
%
% OUTPUT:
%   result         Struct with directions [3 x K], phase angles [rad], speeds
%                  [km/s], delta-v [m/s], flight times [days], perigee altitude
%                  [km], convergence/error arrays, and renderer/provenance
%                  fields. completed does not require every solve to converge.

% Select the isolated historical backend and restore the caller path on exit.
settings = cfg.rpo_to_deep_space.reachability;
overall_timer = tic;
rpo2ds.require_solver(settings.solver_backend);
frozen_directory = fullfile(fileparts(mfilename('fullpath')), 'frozen_reachability');
previous_path = path;
path_cleanup = onCleanup(@() local_restore_path(previous_path));
clear RPO2DS_reachability
addpath(genpath(frozen_directory), '-begin');
expected_solver = fullfile(frozen_directory, 'RPO2DS_reachability.m');
if ~strcmpi(which('RPO2DS_reachability'), expected_solver)
    error('rpo2ds:FrozenReachabilityUnavailable', ...
        'The frozen historical reachability backend could not be selected.');
end
reachability_solver = str2func('RPO2DS_reachability');

% Convert the physical maneuver budget and build the prescribed polar directions.
distance_unit_km = cfg.rpo_to_deep_space.shared.distance_unit_km;
time_unit_seconds = 1 / cfg.rpo_to_deep_space.shared.synodic_angular_rate_rad_s;
velocity_unit_km_s = distance_unit_km / time_unit_seconds;
dv_budget_nd = settings.total_delta_v_budget_m_s / (1000 * velocity_unit_km_s);
angles = linspace(settings.polar_range_rad(1), settings.polar_range_rad(2), ...
    settings.polar_samples).';
directions = zeros(3, numel(angles));
for index = 1:numel(angles)
    directions(:, index) = rotation(settings.initial_direction(:), [0; 1; 0], angles(index));
end
local_validate_preprocessing(preprocessing, settings, directions);

% Solve directions from their archived starts, retaining failures at their grid indices.
result = struct();
result.data = local_empty_data(numel(angles));
timer = tic;
for index = 1:numel(angles)
    row = local_run_point(settings, database, preprocessing.rpo_state_2d, ...
        preprocessing.points(index), directions(:, index), dv_budget_nd, ...
        reachability_solver);
    result.data = local_replace(result.data, row, index);
    if mod(index, max(1, ceil(numel(angles) / 10))) == 0 || index == numel(angles)
        fprintf('  Reachability progress: %d/%d (%.1f s elapsed)\n', ...
            index, numel(angles), toc(timer));
    end
end

% Store physical metrics and the display-only symmetry-reconstruction settings.
result.data.directions = directions;
result.data.polar_angle_rad = angles;
result.data.dv_budget_nd = dv_budget_nd;
result.data.dv_budget_m_s = settings.total_delta_v_budget_m_s;
result.data.renderer = struct('polar_samples', settings.renderer_polar_samples, ...
    'azimuth_samples', settings.renderer_azimuth_samples, ...
    'uses_symmetry_reconstruction', true);
result.configuration = settings;
result.metrics = struct('v_inf_range_km_s', local_range(result.data.v_inf_magnitude_km_s), ...
    'time_of_flight_range_days', local_range(result.data.time_of_flight_days), ...
    'failed_indices', find(~result.data.converged));
result.status = struct('completed', true, 'failed_sections', strings(0, 1));
result.provenance.historical_producer = "Direction_search_optimized_trajs.m / RPO2DS_reachability";
result.provenance.solver_implementation = "frozen pre-refactor fmincon backend";
result.provenance.source_commit = "c9ccc9911c920dcb971746ddbe27fec1e7a6edbb";
result.provenance.preprocessing = preprocessing.provenance;
if ismember('orbitPhase', database.Properties.VariableNames)
    result.provenance.database_role = "phase_aware_generated";
else
    result.provenance.database_role = "legacy_dense_reference";
end
result.provenance.symmetry_reconstructed_surface = true;
result.provenance.elapsed_seconds = toc(overall_timer);
end

function data = local_empty_data(count)
%LOCAL_EMPTY_DATA Allocate ordered directional-sweep arrays and failure diagnostics.

data = struct('directions', zeros(3, count), 'polar_angle_rad', NaN(count, 1), ...
    'v_inf_magnitude_km_s', NaN(count, 1), 'delta_v_m_s', NaN(count, 1), ...
    'time_of_flight_days', NaN(count, 1), 'perigee_altitude_km', NaN(count, 1), ...
    'exitflag', NaN(count, 1), 'converged', false(count, 1), ...
    'objective', NaN(count, 1), 'constraint_residual_inf', NaN(count, 1), ...
    'error_identifier', strings(count, 1), 'error_message', strings(count, 1));
end

function row = local_run_point(settings, database, rpo_state_2d, initial_guess, ...
        direction, dv_budget_nd, solver)
%LOCAL_RUN_POINT Solve one frozen historical direction and record its outcome.
%
% PROTOTYPE:
%   row = local_run_point(settings, database, rpo_state_2d, initial_guess, ...
%       direction, dv_budget_nd, solver)
%
% DESCRIPTION:
%   Pass the archived patch to the selected historical backend. Keep solver
%   status, nonlinear constraint residual, and physical metrics; catch an
%   individual case failure so the remaining directions can still run.
%
% INPUT:
%   settings       Frozen reachability solver/grid settings.
%   database       RPO-Moon table used by the historical solver.
%   rpo_state_2d   [x; y; vx; vy] in rotating DU/VU.
%   initial_guess  Archived patched start for this direction.
%   direction      [3 x 1] Earth-relative direction/seed [km/s].
%   dv_budget_nd   Normalized total impulse budget [VU].
%   solver         Handle to the isolated historical reachability function.
%
% OUTPUT:
%   row            Scalar metric/status/error fields for one polar-grid index.
%                  Speeds use km/s, delta-v m/s, TOF days, altitude km.

row = struct('v_inf_magnitude_km_s', NaN, 'delta_v_m_s', NaN, ...
    'time_of_flight_days', NaN, 'perigee_altitude_km', NaN, 'exitflag', NaN, ...
    'converged', false, 'objective', NaN, 'constraint_residual_inf', NaN, ...
    'error_identifier', "", 'error_message', "");
try
    parameters = [settings.node_count, settings.minimum_earth_altitude_km];
    [v_inf, delta_v, time_of_flight, ~, ~, exitflag, output, altitude] = ...
        feval(solver, direction(:), dv_budget_nd, database, ...
        rpo_state_2d(:), parameters, false, false, initial_guess);
    row.exitflag = exitflag;
    row.objective = -v_inf;
    row.constraint_residual_inf = local_constraint_residual(output);
    % The historical API reports exit flag 0, so feasibility is established
    % from the returned values and the final constraint residual.
    row.converged = isfinite(exitflag) && exitflag >= 0 && ...
        all(isfinite([v_inf, delta_v, time_of_flight, altitude])) && ...
        isfinite(row.constraint_residual_inf) && row.constraint_residual_inf <= 1e-6;
    if row.converged
        row.v_inf_magnitude_km_s = v_inf;
        row.delta_v_m_s = delta_v;
        row.time_of_flight_days = time_of_flight;
        row.perigee_altitude_km = altitude;
    end
catch exception
    row.error_identifier = string(exception.identifier);
    row.error_message = string(exception.message);
end
end

function local_validate_preprocessing(preprocessing, settings, directions)
%LOCAL_VALIDATE_PREPROCESSING Require matching frozen seed/grid/budget preprocessing before the sweep.

% Fail before the long sweep if its frozen preprocessing does not match.

if numel(preprocessing.points) ~= size(directions, 2) || ...
        ~isequal(size(preprocessing.directions), size(directions)) || ...
        max(abs(preprocessing.directions(:) - directions(:))) > 1e-14 || ...
        preprocessing.configuration.polar_samples ~= settings.polar_samples || ...
        preprocessing.configuration.minimum_earth_altitude_km ~= ...
        settings.minimum_earth_altitude_km || ...
        preprocessing.configuration.earth_perigee_samples ~= ...
        settings.earth_perigee_samples || ...
        preprocessing.configuration.rpo_moon_match_tolerance_km_s ~= ...
        settings.rpo_moon_match_tolerance_km_s
    error('rpo2ds:ReachabilityPreprocessingMismatch', ...
        ['The frozen reachability preprocessing does not match the configured ', ...
        '20-direction reference experiment.']);
end
end

function residual = local_constraint_residual(output)
%LOCAL_CONSTRAINT_RESIDUAL Read the final historical fmincon constraint residual, or NaN.

% Read the final historical fmincon constraint residual.

residual = NaN;
if iscell(output) && ~isempty(output)
    final_output = output{end};
    if isstruct(final_output) && isfield(final_output, 'constrviolation')
        residual = final_output.constrviolation;
    end
end
end

function local_restore_path(previous_path)
%LOCAL_RESTORE_PATH Restore the caller path and clear cached historical solver functions.

% Remove the isolated historical backend after this section completes.

path(previous_path);
clear RPO2DS_reachability constr_CR3BP_reachability objfun_CR3BP_reachability ...
    earth_flyby_for_optimization get_dr1ddot_dtX get_dr2ddot_dtX ...
    get_r1ddot get_r2ddot propagate distance_from_point obj_zero ...
    RTH2ECI_moon in2rot rot2in event_xaxiscrossing astroConstants ode_2bp ...
    rhs_BCR4BP_2D rhs_BCR4BP_2D_STM rhs_CR3BP_2D rhs_CR3BP_2D_STM ...
    rhs_CR3BP_3D rhs_CR3BP_3D_STM
end

function data = local_replace(data, row, index)
%LOCAL_REPLACE Store one direction result without changing its polar-grid index.

data.v_inf_magnitude_km_s(index) = row.v_inf_magnitude_km_s;
data.delta_v_m_s(index) = row.delta_v_m_s;
data.time_of_flight_days(index) = row.time_of_flight_days;
data.perigee_altitude_km(index) = row.perigee_altitude_km;
data.exitflag(index) = row.exitflag;
data.converged(index) = row.converged;
data.objective(index) = row.objective;
data.constraint_residual_inf(index) = row.constraint_residual_inf;
data.error_identifier(index) = row.error_identifier;
data.error_message(index) = row.error_message;
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
