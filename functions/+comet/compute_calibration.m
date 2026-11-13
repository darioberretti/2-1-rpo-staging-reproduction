function result = compute_calibration(settings,database,seed)
%COMPUTE_CALIBRATION Calibrate the physical RPO budget and solve a planar RTH capability boundary.
%
% PROTOTYPE:
%   result = comet.compute_calibration(settings, database, seed)
%
% DESCRIPTION:
%   Retain the historical endpoint, bisection and interpolation
%   sequence from RPO_planar_sweep_for_sanchez.m. Failed intermediate
%   trial costs remain in that search, but the closest-cost anchor must
%   pass exitflag, explicit constraint/bound and physical-budget checks.
%   Use its regularized budget for a serial directional sweep, requiring
%   every stage-2 solve to pass before returning a complete boundary.
%   Monte Carlo uses its separately configured published scalar limit.
%
% INPUT:
%   settings  cfg.comet_case_study.calibration from reproduction_config:
%               first (direction [deg], budget/tolerance [m/s], excess-speed
%               bracket [km/s], iteration/continuation controls), sweep
%               (angle grid [deg], seed/solver controls), params/opts
%               (RPO2DS_transfer shooting and optimizer settings), and log.
%   database  Phase-aware RPO-Moon table from rpo2ds.load_database;
%               its schema/normalization is documented in docs/DATA.md.
%   seed      [4 elements] Planar rotating EM state [x,y,vx,vy], DU/VU,
%               from data/inputs/reference_rpo_seed.mat.
%
% OUTPUT:
%   result    Struct with settings, calibration (history table with
%               vinf_mag_km_s [km/s], dv_m_s [m/s], regularized_budget_nd,
%               tof_days [days], exitflag, method and constraint_violation;
%               best_row, direction_rth and compact trial solutions),
%               sweep (angle/direction, speed/cost/time and solver table),
%               anchor, solutions (compact_solution schema), budget_nd
%               (internal normalized sum of squared impulse variables)
%               and elapsed_seconds [s]. Directions use Earth-Sun RTH.
%               Source/seed hashes and common provenance are attached
%               by run_06_comet_case_study after this function returns.

rpo2ds.require_solver(settings.first.solver_backend);
rpo2ds.require_solver(settings.sweep.solver_backend);

% Calibrate the first direction and reject an unacceptable final anchor.
timer = tic;
direction = local_planar_rth_direction(settings.first.direction_angle_deg);
[history,anchor] = local_calibrate_first_direction(direction,database,seed,settings);
if anchor.exitflag<=0 || anchor.calibration_constraint_violation>settings.opts.constraint_tolerance ...
        || abs(anchor.dv_opt-settings.first.requested_budget_m_s)>settings.first.budget_tolerance_m_s
    error('comet:Calibration','RPO calibration did not converge within the requested maneuver-budget tolerance. The directional sweep was not started.');
end

% Use the validated anchor budget for the serial directional continuation.
budget = local_regularized_budget_nd(anchor);
[sweep,solutions] = local_run_planar_reachability_sweep(anchor,budget,database,seed,settings);

% Retain complete trial diagnostics alongside compact accepted sweep solutions.
result = struct('settings',settings,'calibration',history,'sweep',sweep, ...
    'anchor',compact_solution(anchor),'solutions',{solutions}, ...
    'budget_nd',budget,'elapsed_seconds',toc(timer));
fprintf('RPO calibration anchor %.3f m/s; boundary %.6f--%.6f km/s; %.1f s\n', ...
    anchor.dv_opt,min(sweep.reachable_vinf_mag_km_s),max(sweep.reachable_vinf_mag_km_s),toc(timer));
end

function opts = local_apply_solver_controls(opts, solver_settings)
%LOCAL_APPLY_SOLVER_CONTROLS Copy the selected backend, stage-1 toggle and shared algorithm into opts.

opts.solver_backend = solver_settings.solver_backend;
opts.run_stage1 = solver_settings.run_stage1;
opts.stage1_algorithm = solver_settings.fmincon_algorithm;
opts.stage2_algorithm = solver_settings.fmincon_algorithm;
end

function [budget_continuation, best_result] = local_calibrate_first_direction(direction_rth, RPOMoonDatabase, xx0_RPO, settings)
%LOCAL_CALIBRATE_FIRST_DIRECTION Tune the first fixed-v-infinity direction with the historical budget search.
%
% PROTOTYPE:
%   [budget_continuation, best_result] = local_calibrate_first_direction( ...
%       direction_rth, RPOMoonDatabase, xx0_RPO, settings)
%
% DESCRIPTION:
%   Solve bracket endpoints, bisect while needed, then interpolate.
%   Sort active rows by target excess-speed magnitude after each solve.
%   Select the returned physical cost closest to the target; final
%   anchor acceptance is enforced by compute_calibration, not here.
%
% INPUT:
%   direction_rth   [3 x 1] Unit Earth-Sun RTH direction [-].
%   RPOMoonDatabase Phase-aware encounter table for RPO2DS_transfer.
%   xx0_RPO         [4 elements] Planar EM rotating seed [x,y,vx,vy], DU/VU.
%   settings        Calibration settings described in compute_calibration.
%
% OUTPUT:
%   budget_continuation  Struct with direction, requested_budget_m_s,
%                         full history table, best_row and compact solutions.
%   best_result          Raw RPO2DS_transfer result with the closest cost,
%                         including calibration_constraint_violation.

params = settings.params;
opts = settings.opts;
opts.graph = false; % Keep the calibration loop focused on numeric outputs.
opts.setup_only = false; % Force actual solves during the calibration loop.
opts = local_apply_solver_controls(opts, settings.first); % Apply the first-direction solver controls.
target_budget = settings.first.requested_budget_m_s; % Requested physical RPO budget [m/s]
target_tolerance = settings.first.budget_tolerance_m_s; % Accepted budget error [m/s]
bracket_width_for_interp = settings.first.interpolation_switch_width_km_s; % Width that enables interpolation [km/s]

% Preserve the original endpoint/bisection/interpolation sequence.
initial_mags = settings.first.vinf_mag_bracket_km_s(:);
max_history_rows = numel(initial_mags) + settings.first.bisection_iterations + settings.first.interpolation_iterations;
mag_history = NaN(max_history_rows, 1);
result_history = cell(max_history_rows, 1);
dv_history = NaN(max_history_rows, 1);
budget_history = NaN(max_history_rows, 1);
tof_history = NaN(max_history_rows, 1);
exit_history = NaN(max_history_rows, 1);
method_history = strings(max_history_rows, 1);
n_history = 0;
base_result = [];
use_calibration_continuation = logical(settings.first.use_continuation); % Toggle warm-start reuse for calibration solves.

% Solve the two bracket endpoints independently.
for i_endpoint = 1:numel(initial_mags)
    n_history = n_history + 1;
    mag_history(n_history) = initial_mags(i_endpoint);
    method_history(n_history) = "bracket";
    [result_history{n_history}, dv_history(n_history), budget_history(n_history), ...
        tof_history(n_history), exit_history(n_history), base_result] = local_solve_calibration_point( ...
        mag_history(n_history), direction_rth, RPOMoonDatabase, xx0_RPO, params, opts, ...
        base_result, use_calibration_continuation, settings, n_history, max_history_rows, method_history(n_history));
end

% Keep the history sorted so neighboring rows define candidate brackets.
[mag_history, result_history, dv_history, budget_history, tof_history, exit_history, method_history] = ...
    local_sort_calibration_history(mag_history, result_history, dv_history, budget_history, tof_history, exit_history, method_history, n_history);

% Use bisection until the bracket is tight enough or the target is reached.
for i_bisect = 1:settings.first.bisection_iterations
    [bracket_idx, closest_idx] = local_find_budget_bracket(dv_history(1:n_history), target_budget);
    if abs(dv_history(closest_idx) - target_budget) <= target_tolerance
        break
    end
    if isempty(bracket_idx)
        break
    end
    mag_a = mag_history(bracket_idx);
    mag_b = mag_history(bracket_idx + 1);
    if abs(mag_b - mag_a) <= bracket_width_for_interp
        break
    end
    n_history = n_history + 1;
    mag_history(n_history) = 0.5 * (mag_a + mag_b);
    method_history(n_history) = "bisection";
    [result_history{n_history}, dv_history(n_history), budget_history(n_history), ...
        tof_history(n_history), exit_history(n_history), base_result] = local_solve_calibration_point( ...
        mag_history(n_history), direction_rth, RPOMoonDatabase, xx0_RPO, params, opts, ...
        base_result, use_calibration_continuation, settings, n_history, max_history_rows, method_history(n_history));
    [mag_history, result_history, dv_history, budget_history, tof_history, exit_history, method_history] = ...
        local_sort_calibration_history(mag_history, result_history, dv_history, budget_history, tof_history, exit_history, method_history, n_history);
end

% Use interpolation after bisection has made a local bracket.
for i_interp = 1:settings.first.interpolation_iterations
    [bracket_idx, closest_idx] = local_find_budget_bracket(dv_history(1:n_history), target_budget);
    if abs(dv_history(closest_idx) - target_budget) <= target_tolerance
        break
    end
    if isempty(bracket_idx)
        break
    end
    mag_a = mag_history(bracket_idx);
    mag_b = mag_history(bracket_idx + 1);
    dv_a = dv_history(bracket_idx);
    dv_b = dv_history(bracket_idx + 1);
    if abs(dv_b - dv_a) <= eps
        break
    end
    n_history = n_history + 1;
    mag_history(n_history) = mag_a + (target_budget - dv_a) * (mag_b - mag_a) / (dv_b - dv_a);
    method_history(n_history) = "interpolation";
    [result_history{n_history}, dv_history(n_history), budget_history(n_history), ...
        tof_history(n_history), exit_history(n_history), base_result] = local_solve_calibration_point( ...
        mag_history(n_history), direction_rth, RPOMoonDatabase, xx0_RPO, params, opts, ...
        base_result, use_calibration_continuation, settings, n_history, max_history_rows, method_history(n_history));
    [mag_history, result_history, dv_history, budget_history, tof_history, exit_history, method_history] = ...
        local_sort_calibration_history(mag_history, result_history, dv_history, budget_history, tof_history, exit_history, method_history, n_history);
end

% Select the closest physical-budget solve as the reachability budget anchor.
mag_history = mag_history(1:n_history);
result_history = result_history(1:n_history);
dv_history = dv_history(1:n_history);
budget_history = budget_history(1:n_history);
tof_history = tof_history(1:n_history);
exit_history = exit_history(1:n_history);
method_history = method_history(1:n_history);
[~, best_idx] = min(abs(dv_history - target_budget));
best_result = result_history{best_idx};
best_row = table(mag_history(best_idx), dv_history(best_idx), budget_history(best_idx), ...
    tof_history(best_idx), exit_history(best_idx), ...
    'VariableNames', {'vinf_mag_km_s', 'dv_m_s', 'regularized_budget_nd', 'tof_days', 'exitflag'});

% Store the complete continuation history for plotting and later inspection.
budget_continuation = struct();
budget_continuation.direction_rth = direction_rth;
budget_continuation.requested_budget_m_s = target_budget;
budget_continuation.history = table(mag_history, dv_history, budget_history, tof_history, exit_history, method_history, ...
    'VariableNames', {'vinf_mag_km_s', 'dv_m_s', 'regularized_budget_nd', 'tof_days', 'exitflag', 'method'});
budget_continuation.history.constraint_violation = cellfun( ...
    @(raw)raw.calibration_constraint_violation,result_history);
budget_continuation.best_row = best_row;
budget_continuation.solutions = cellfun(@compact_solution,result_history,'UniformOutput',false);
end

function [result_i, dv_i, budget_i, tof_i, exit_i, base_result] = local_solve_calibration_point( ...
    vinf_mag_km_s, direction_rth, RPOMoonDatabase, xx0_RPO, params, opts, ...
    base_result, use_continuation, settings, solve_idx, max_solves, solve_method)
%LOCAL_SOLVE_CALIBRATION_POINT Solve one fixed-speed calibration point and retain its measured feasibility.
%
% PROTOTYPE:
%   [result_i, dv_i, budget_i, tof_i, exit_i, base_result] = ...
%       local_solve_calibration_point(vinf_mag_km_s, direction_rth, ...
%       RPOMoonDatabase, xx0_RPO, params, opts, base_result, ...
%       use_continuation, settings, solve_idx, max_solves, solve_method)
%
% DESCRIPTION:
%   Optionally warm-start from the previous raw solve, then evaluate
%   nonlinear constraints and decision bounds explicitly. Warn on a failed
%   trial while retaining its returned cost for the historical search.
%
% INPUT:
%   vinf_mag_km_s   Target excess-speed magnitude [km/s].
%   direction_rth   [3 x 1] Unit Earth-Sun RTH direction [-].
%   RPOMoonDatabase, xx0_RPO  Encounter table and normalized planar seed.
%   params, opts    RPO2DS_transfer parameter and optimizer structs.
%   base_result     Previous raw result, or [] before the first solve.
%   use_continuation Logical enabling reuse when base_result is nonempty.
%   settings        Calibration settings, including tolerance and log.
%   solve_idx, max_solves, solve_method  Progress counter, cap and label.
%
% OUTPUT:
%   result_i   Raw result plus maximum normalized constraint/bound violation.
%   dv_i       Physical optimized impulse cost [m/s].
%   budget_i   Sum of squared regularized impulse variables, normalized.
%   tof_i      Flight time [days].
%   exit_i     Returned optimizer status.
%   base_result  Current result for a possible later warm start.

params.vinf_target_rth_km_s = vinf_mag_km_s * direction_rth;
opts_i = opts;
has_continuation = use_continuation && ~isempty(base_result);
if has_continuation
    opts_i.continuation_base_result = base_result;
end
local_progress_log(settings.log, 'SANCHEZ', 'calibration %s solve %d/%d: |v_inf| %.6f km/s, continuation %d', solve_method, solve_idx, max_solves, vinf_mag_km_s, has_continuation);
result_i = RPO2DS_transfer('opt', params.vinf_target_rth_km_s, RPOMoonDatabase, xx0_RPO, params, opts_i);
[c,ceq] = constr_CR3BP_regularized_phase(result_i.y_opt,result_i.setup.problem_data);
result_i.calibration_constraint_violation = max([0;c(:);abs(ceq(:)); ...
    result_i.setup.lb(:)-result_i.y_opt(:);result_i.y_opt(:)-result_i.setup.ub(:)]);
if result_i.exitflag<=0 || result_i.calibration_constraint_violation>settings.opts.constraint_tolerance
    % The original budget search uses returned trial costs even when a trial
    % fails. Keep its search sequence, but require a valid final anchor above.
    warning('comet:CalibrationPoint','Calibration trial at %.6f km/s did not converge (exit %g; constraint/bound error %.3e). Its returned cost is retained for the historical bracket search; this is not a validated trajectory.', ...
        vinf_mag_km_s,result_i.exitflag,result_i.calibration_constraint_violation);
end
dv_i = result_i.dv_opt;
budget_i = local_regularized_budget_nd(result_i);
tof_i = result_i.tof_opt;
exit_i = result_i.exitflag;
base_result = result_i;
local_progress_log(settings.log, 'SANCHEZ', 'calibration %s solve %d complete: dv %.6f m/s, tof %.6f days, exit %g', solve_method, solve_idx, dv_i, tof_i, exit_i);
end

function [bracket_idx, closest_idx] = local_find_budget_bracket(dv_history, target_budget)
%LOCAL_FIND_BUDGET_BRACKET Find the first adjacent cost pair straddling target_budget [m/s].
% Return its first index (or []), plus the closest-cost row index;
% dv_history must follow the magnitude-sorted calibration history.

diff_history = dv_history - target_budget;
[~, closest_idx] = min(abs(diff_history));
bracket_idx = find(diff_history(1:end-1) .* diff_history(2:end) <= 0, 1, 'first');
end

function [mag_history, result_history, dv_history, budget_history, tof_history, exit_history, method_history] = ...
    local_sort_calibration_history(mag_history, result_history, dv_history, budget_history, tof_history, exit_history, method_history, n_history)
%LOCAL_SORT_CALIBRATION_HISTORY Sort active rows by target speed and apply the same order to every history.
% The first n_history entries are active; preallocated trailing rows remain.
% Magnitude uses km/s, cost m/s, TOF days; results/status/method follow rows.

[~, sort_idx] = sort(mag_history(1:n_history));
mag_history(1:n_history) = mag_history(sort_idx);
result_history(1:n_history) = result_history(sort_idx);
dv_history(1:n_history) = dv_history(sort_idx);
budget_history(1:n_history) = budget_history(sort_idx);
tof_history(1:n_history) = tof_history(sort_idx);
exit_history(1:n_history) = exit_history(sort_idx);
method_history(1:n_history) = method_history(sort_idx);
end

function [sweep_table, sweep_results] = local_run_planar_reachability_sweep(first_solution, reach_budget_nd, RPOMoonDatabase, xx0_RPO, settings)
%LOCAL_RUN_PLANAR_REACHABILITY_SWEEP Solve every planar Earth-Sun RTH direction at the calibrated internal budget.
%
% PROTOTYPE:
%   [sweep_table, sweep_results] = local_run_planar_reachability_sweep( ...
%       first_solution, reach_budget_nd, RPOMoonDatabase, xx0_RPO, settings)
%
% DESCRIPTION:
%   Promote the anchor decision vector to reachability mode, adding
%   a free target-speed magnitude. With continuation enabled, rotate the
%   free Sun-phase seed to preserve the inertial target as RTH direction
%   changes. Abort on any failed stage-2 exitflag/constraint check.
%
% INPUT:
%   first_solution  Validated raw anchor result from RPO2DS_transfer.
%   reach_budget_nd Normalized sum-of-squared-impulse budget from the anchor.
%   RPOMoonDatabase Phase-aware RPO-Moon encounter table.
%   xx0_RPO         [4 elements] Planar EM rotating seed [x,y,vx,vy], DU/VU.
%   settings        Calibration settings with sweep.angle_grid_deg,
%                     solve_mode (continuation enables vector reuse),
%                     use_previous_magnitude_seed, params and opts.
%
% OUTPUT:
%   sweep_table     One row per angle: angle_deg [deg], unit R/T/H
%                     components, seed/reachable_vinf_mag_km_s [km/s],
%                     solve mode, dv_m_s [m/s], regularized_budget_nd,
%                     tof_days [days], theta_sun0_rad [rad] and exitflag.
%   sweep_results   Column cell array of compact_solution records.

angle_grid = settings.sweep.angle_grid_deg(:); % Sweep support [deg]
n_angles = numel(angle_grid); % Number of solved direction samples
use_sweep_continuation = strcmpi(settings.sweep.solve_mode, 'continuation'); % True reuses the previous optimized decision vector
params = settings.params;
opts = settings.opts;
opts.graph = false; % Disable per-angle trajectory plots during the sweep.
opts.setup_only = false; % Force actual reachability solves.
opts.dv_budget = reach_budget_nd; % Internal reachability budget calibrated from the first direction.
opts = local_apply_solver_controls(opts, settings.sweep); % Apply the continuation sweep solver controls.

% Preallocate scalar sweep diagnostics.
vInfMag_km_s = NaN(n_angles, 1);
dv_m_s = NaN(n_angles, 1);
tof_days = NaN(n_angles, 1);
exitflag = NaN(n_angles, 1);
theta_sun0_rad = NaN(n_angles, 1);
regularized_budget_nd = NaN(n_angles, 1);
direction_R = NaN(n_angles, 1);
direction_T = NaN(n_angles, 1);
direction_H = NaN(n_angles, 1);
vinf_seed_mag_km_s = NaN(n_angles, 1);
sweep_solve_mode = strings(n_angles, 1);
sweep_results = cell(n_angles, 1);

% Start the first direction from the validated calibration trajectory.
previous_result = first_solution;
previous_result.y_opt = [first_solution.y_opt(:);norm(first_solution.vinf_target_rth)];
previous_result.setup.problem_data.mode = 'reach';
previous_result.setup.problem_data.vinf_mag_idx = numel(previous_result.y_opt);
previous_result.setup.problem_data.nvars = numel(previous_result.y_opt);
previous_result.setup.problem_data.dv_budget = reach_budget_nd;
previous_result.setup.problem_data.vinf_dir = first_solution.vinf_target_rth(:)/norm(first_solution.vinf_target_rth);
previous_result.setup.lb = [first_solution.setup.lb(:);0];
previous_result.setup.ub = [first_solution.setup.ub(:);Inf];
previous_result.meta.vinf_mag_idx = numel(previous_result.y_opt);
previous_result.meta.decision_vector_length = numel(previous_result.y_opt);
previous_result.meta.vInfMag0 = norm(first_solution.vinf_target_rth);
previous_result.meta.dv_budget = reach_budget_nd;
previous_result.meta.vinf_dir = previous_result.setup.problem_data.vinf_dir;
previous_mag = norm(first_solution.vinf_target_rth);

% March around the planar circle with the selected sweep solve mode.
for i_angle = 1:n_angles
    direction_rth = local_planar_rth_direction(angle_grid(i_angle));
    direction_R(i_angle) = direction_rth(1);
    direction_T(i_angle) = direction_rth(2);
    direction_H(i_angle) = direction_rth(3);
    vinf_seed_mag_km_s(i_angle) = previous_mag;
    sweep_solve_mode(i_angle) = string(settings.sweep.solve_mode);
    params.vinf_target_rth_km_s = previous_mag * direction_rth;
    opts_i = opts;
    has_continuation = use_sweep_continuation && ~isempty(previous_result);
    if has_continuation
        opts_i.continuation_base_result = previous_result;
        % Free Sun phase can absorb the RTH direction change while retaining
        % the same physical ECI target and a feasible trajectory seed.
        previous_direction = previous_result.vinf_target_rth(:);
        angle_change = atan2(direction_rth(2),direction_rth(1)) ...
            - atan2(previous_direction(2),previous_direction(1));
        opts_i.y0_override = previous_result.y_opt;
        sun_phase_index = previous_result.meta.theta_sun0_idx;
        opts_i.y0_override(sun_phase_index) = mod( ...
            opts_i.y0_override(sun_phase_index)-angle_change,2*pi);
    end
    local_progress_log(settings.log, 'SANCHEZ', 'sweep solve %d/%d start: angle %.3f deg, seed |v_inf| %.6f km/s, mode %s, continuation %d', i_angle, n_angles, angle_grid(i_angle), previous_mag, settings.sweep.solve_mode, has_continuation);
    result_i = RPO2DS_transfer('reach', params.vinf_target_rth_km_s, RPOMoonDatabase, xx0_RPO, params, opts_i);
    if result_i.solver.exitflag_stage2<=0 || result_i.solver.stage2.constrviolation>settings.opts.constraint_tolerance
        error('comet:DirectionalSolve','RPO boundary solve at %.0f degrees did not converge (exit %g; constraint error %.3e). No complete boundary was saved.', ...
            angle_grid(i_angle),result_i.solver.exitflag_stage2,result_i.solver.stage2.constrviolation);
    end
    result_i.exitflag = result_i.solver.exitflag_stage2; % The wrapper's reach flag is a legacy placeholder.
    result_i = strip_previous(result_i);
    sweep_results{i_angle} = compact_solution(result_i);
    vInfMag_km_s(i_angle) = result_i.vInfMag;
    dv_m_s(i_angle) = result_i.dv_opt;
    tof_days(i_angle) = result_i.tof_opt;
    exitflag(i_angle) = result_i.exitflag;
    regularized_budget_nd(i_angle) = local_regularized_budget_nd(result_i);
    theta_sun0_rad(i_angle) = local_read_theta_sun0(result_i);
    if settings.sweep.use_previous_magnitude_seed && isfinite(result_i.vInfMag)
        previous_mag = result_i.vInfMag;
    end
    previous_result = result_i;
    fprintf('Sweep %4d/%4d | angle %8.3f deg | |v_inf| %.6f km/s | dv %.3f m/s\n', ...
        i_angle, n_angles, angle_grid(i_angle), vInfMag_km_s(i_angle), dv_m_s(i_angle));
    drawnow('limitrate');
    local_progress_log(settings.log, 'SANCHEZ', 'sweep solve %d complete: |v_inf| %.6f km/s, dv %.6f m/s, tof %.6f days, exit %g', i_angle, vInfMag_km_s(i_angle), dv_m_s(i_angle), tof_days(i_angle), exitflag(i_angle));
end

% Package the scalar boundary summary for export and capability inspection.
sweep_table = table(angle_grid, direction_R, direction_T, direction_H, vinf_seed_mag_km_s, ...
    sweep_solve_mode, vInfMag_km_s, dv_m_s, regularized_budget_nd, tof_days, theta_sun0_rad, exitflag, ...
    'VariableNames', {'angle_deg', 'direction_R', 'direction_T', 'direction_H', ...
    'vinf_seed_mag_km_s', 'sweep_solve_mode', 'reachable_vinf_mag_km_s', 'dv_m_s', ...
    'regularized_budget_nd', 'tof_days', 'theta_sun0_rad', 'exitflag'});
end

function direction_rth = local_planar_rth_direction(angle_deg)
%LOCAL_PLANAR_RTH_DIRECTION Return [cos(angle);sin(angle);0], a [3 x 1] Earth-Sun RTH unit vector.
% angle_deg is measured from +R toward +T in degrees.

direction_rth = [cosd(angle_deg); sind(angle_deg); 0];
direction_rth = direction_rth / norm(direction_rth);
end

function regularized_budget_nd = local_regularized_budget_nd(result)
%LOCAL_REGULARIZED_BUDGET_ND Sum squared impulse variables from the first N nine-entry node blocks.
% Each block is [x,y,z,vx,vy,vz,u,w,s]; return the normalized reach budget.

n_nodes = result.nodes.n_nodes;
node_blocks = reshape(result.y_opt(1:9 * n_nodes), 9, n_nodes);
regularized_budget_nd = sum(sum(node_blocks(7:9, :).^2));
end

function theta_sun0 = local_read_theta_sun0(result)
%LOCAL_READ_THETA_SUN0 Read optimized Sun phase [rad] via meta.theta_sun0_idx, or return NaN.

theta_sun0 = NaN;
if isfield(result, 'meta') && isfield(result.meta, 'theta_sun0_idx') && ~isempty(result.meta.theta_sun0_idx)
    theta_sun0 = result.y_opt(result.meta.theta_sun0_idx);
end
end

function local_progress_log(~,component,message,varargin)
%LOCAL_PROGRESS_LOG Print timestamped sprintf progress and process pending graphics events.
% The first compatibility argument is ignored; logging is unconditional.

fprintf('[%s] %s: %s\n',char(datetime('now','Format','HH:mm:ss')),component,sprintf(message,varargin{:}));
drawnow('limitrate');
end

function raw = strip_previous(raw)
%STRIP_PREVIOUS Remove the nested continuation_base_result from raw.inputs.opts, if present.

if isfield(raw.inputs.opts,'continuation_base_result')
    raw.inputs.opts = rmfield(raw.inputs.opts,'continuation_base_result');
end
end

function solution = compact_solution(raw)
%COMPACT_SOLUTION Keep trajectory/solver evidence without recursive optimizer inputs.
%
% PROTOTYPE:
%   solution = compact_solution(raw)
%
% DESCRIPTION:
%   Keep one representative trajectory entry because the wrapper
%   expands cell-valued STMs into repeated struct entries carrying the
%   same node arrays. Retain decision bounds and problem data for checks.
%
% INPUT:
%   raw       RPO2DS_transfer result with optimizer outputs and, for
%               a CasADi trial, calibration_constraint_violation.
%
% OUTPUT:
%   solution  Struct with y_opt, dv_m_s [m/s], tof_days [days], exitflag,
%               constants, meta, phase_final [TU], nodes, impulses, dv,
%               trajectory, solver_output, vinf_target_rth_kms [3 x 1, km/s],
%               problem_data and bounds.lower/upper. trajectory.t_nodes
%               uses EM TU; xx_nodes [N x 6] uses rotating EM DU/VU in
%               [x,y,z,vx,vy,vz] order; U_nodes [N x 3] stores regularized
%               impulses. Solver output records backend, iterations,
%               message, constraint_violation and first_order_optimality.
%               Full raw interface: RPO2DS_transfer and its local build_result_struct.

solution = struct('y_opt',raw.y_opt,'dv_m_s',raw.dv_opt, ...
    'tof_days',raw.tof_opt,'exitflag',raw.exitflag,'constants',raw.constants, ...
    'meta',raw.meta,'phase_final',raw.phase_final);
solution.nodes = raw.nodes;
solution.impulses = raw.impulses;
% Full postprocessing expands the cell-valued STMs into a struct array;
% every entry carries the same node arrays. Only these arrays are retained.
trajectory = raw.trajectory(1);
solution.trajectory = struct('t_nodes',trajectory.t_nodes, ...
    'xx_nodes',trajectory.xx_nodes,'U_nodes',trajectory.U_nodes);
solution.dv = raw.dv;
output = raw.solver.stage2;
if isfield(output,'return_status')
    solution.solver_output = struct('backend',"casadi_ipopt", ...
        'iterations',output.stats.iter_count,'message',output.return_status, ...
        'constraint_violation',raw.calibration_constraint_violation, ...
        'first_order_optimality',output.stats.iterations.inf_du(end));
else
    solution.solver_output = struct('backend',"fmincon", ...
        'iterations',output.iterations,'message',output.message, ...
        'constraint_violation',output.constrviolation,'first_order_optimality',output.firstorderopt);
end
solution.vinf_target_rth_kms = raw.vinf_target_rth;
solution.problem_data = raw.setup.problem_data;
solution.bounds = struct('lower',raw.setup.lb,'upper',raw.setup.ub);
end
