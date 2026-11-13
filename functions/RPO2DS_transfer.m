function result = RPO2DS_transfer(mode, vInfPlusEarth, RPOMoonDatabase, xx0_RPO, params, opts)
%RPO2DS_TRANSFER Unified transfer optimization/reachability solver.
%
% PROTOTYPE:
%   result = RPO2DS_transfer(mode, vInfPlusEarth, RPOMoonDatabase, xx0_RPO, params, opts)
%
% DESCRIPTION:
%   Build a patched RPO-Moon-Earth seed, reconstruct the lunar flyby, place
%   multiple-shooting nodes, and run the selected feasibility/objective solve.
%   Setup-only and continuation paths use the same result structures. Run02
%   uses the frozen CasADi/IPOPT validation profile; historical reachability
%   is handled by its separate backend under +rpo2ds/frozen_reachability.
%
%   Rotating states are normalized [x; y; z; vx; vy; vz] in DU/VU, with
%   absolute epochs and durations in TU. Physical patch tables use km, km/s,
%   seconds/days as indicated by their fields. RPO phase is a time offset
%   in TU, while Sun/Moon geometry angles are in radians.
%
% INPUT:
%   mode           : 'opt'   -> minimize total DV (regularized)
%                    'reach' -> maximize vInf magnitude under DV budget
%   vInfPlusEarth  : 3x1 vector (km/s) for target vInf direction/magnitude
%   RPOMoonDatabase: table used by RPO_Moon_Earth_guess
%   xx0_RPO        : initial RPO state (4x1, planar)
%   params         : struct (all optional, defaults apply)
%       .n_nodes               (default 40)
%       .hStart                (default 1000 km)
%       .N                     (default 101)
%       .sweep_rad             (default 2*pi)
%       .nSol                  (default 1)
%       .feasibleLimit         (default 1e-2)
%       .maxTOF                (default 100 days)
%       .maxTransferTOF        transfer-only cap in days used by the
%                               operational modes (default maxTOF)
%       .nMoonOrbitSamples     (default 100 full-circle Moon-phase samples)
%       .time_mode             (default 'optimization_only')
%       .maxTOF_relax_factor   multiplier applied to one selected over-cap
%                               max_time seed before enforcing the NLP TOF
%                               limit (default 1.1)
%       .time_variable_mode    (default 'node_times')
%                             : 'node_times' packs absolute node times,
%                               'propagation_times' packs dt between nodes
%       .geometry_mode         (default 'free_geometry')
%                             : 'free_geometry' leaves node_times unchanged and
%                               adds t(1) to propagation_times
%                               'prescribed_geometry' uses one prescribed
%                               initial geometry epoch when selecting the
%                               Moon-Earth / RPO-Moon patch candidate and
%                               pins t(1) to that same epoch in the NLP
%       .initial_geometry_time prescribed initial geometry epoch in
%                               nondimensional TU used only when
%                               geometry_mode is 'prescribed_geometry'
%       .t_node1_fixed         inactive legacy field kept for compatibility
%       .t_nodeN_fixed         inactive legacy field kept for compatibility
%       .dense_nodes_multiplier(default 0.15)
%       .node_count_mode       (default 'fixed')
%       .nodes_per_day         (used when node_count_mode is 'adaptive' or 'adaptive_plus')
%       .n_nodes_min           (adaptive/adaptive_plus lower clamp)
%       .n_nodes_max           (adaptive/adaptive_plus upper clamp)
%       .adaptive_plus_divisions     (default 6)
%       .adaptive_plus_upstream_nodes(default 2)
%       .phase_mode            (default 'seed_phase')
%       .phase_value           (used only when phase_mode is 'other_phase')
%       .impulse_mode          (default 'every_node')
%                             : 'every_node' allows burns at every node,
%                               'endpoints_only' pins interior node burns to zero
%   opts           : struct (all optional, defaults apply)
%       .graph            struct of plot flags:
%                          flyby, initial_guess, nodes,
%                          optimized_rot, optimized_inertial
%       .dv_budget        required for 'reach' (nondimensional)
%       .gradient_check   run checkGradients and return early (default 0)
%       .jacobian_backend 'dense' (default) or 'sparse'
%       .solver_display   fmincon display level (default 'iter')
%       .solver_backend   'fmincon' (default) or 'casadi_ipopt'
%       .casadi_path      CasADi MATLAB interface folder for IPOPT backend
%       .print_logs       print compact setup and solve logs (default true)
%       .print_patch_logs print max_time patch-search diagnostics (default true)
%       .patch_log_candidate_rows maximum patch candidate rows printed (default Inf)
%       .setup_only       stop after building the initial guess and NLP data
%
% OUTPUT: result struct
%   .dv_opt, .tof_opt, .y_opt, .J_opt, .exitflag, .output, .h
%   .vInfMag (reach mode only)
%   .time_variable_mode active NLP time layout used by setup/solve/postprocess
%   .FullSolution and rich sub-structs:
%       inputs, constants, initial_guess, nodes, impulses,
%       solver, trajectory, dv, plots, meta
%   .trajectory.t_nodes always stores absolute node times, even when the
%       decision vector uses propagation durations
%   .setup carries y0/y0_seed, lb/ub, and problem_data for the active NLP.
%       Node blocks occupy the first 9*n variables: [X(6); U(3)] per node.
%       U regularizes the impulse, with norm(delta-v) = u^2 + w^2 + s^2.
%       problem_data records time/phase/target indices and model constants;
%       .meta records decoding and geometry conventions for postprocessing.
%   .constants gives DU [km], TU [s], TU_days [days], VU [km/s], and mu_EM.
%   .trajectory rotating coordinates use DU/VU; inertial positions use km.
%       .dv.DVvec_opt_in uses m/s; dv_opt is in m/s and tof_opt in days.
%

if nargin < 5 || isempty(params)
    params = struct();
end
if nargin < 6 || isempty(opts)
    opts = struct();
end

[params, opts] = apply_defaults(params, opts);
local_progress_log(opts.progress_logs, 'RPO2DS', 'defaults applied; requested mode %s, backend %s', char(string(mode)), char(opts.solver_backend));

mode = lower(string(mode));
if opts.solver_backend == "casadi_ipopt"
    local_progress_log(opts.progress_logs, 'RPO2DS', 'applying CasADi/IPOPT profile');
    [mode, params, opts] = local_apply_casadi_ipopt_profile_v1(mode, params, opts);
    local_progress_log(opts.progress_logs, 'RPO2DS', 'CasADi/IPOPT profile applied; active mode %s', char(mode));
end
params.operational_mode = lower(string(params.operational_mode));
params.phase_mode = lower(string(params.phase_mode));
params.time_mode = lower(string(params.time_mode));
    params.time_variable_mode = lower(string(params.time_variable_mode));
    params.geometry_mode = lower(string(params.geometry_mode));
    params.target_frame = lower(string(params.target_frame));
    params.theta_sun0_mode = lower(string(params.theta_sun0_mode));
    params.final_time_mode = lower(string(params.final_time_mode));
    if params.operational_mode == "fixed_arrival"
        params.final_time_mode = "specified";
        params.final_time_seed = params.final_epoch_time;
    else
        params.final_time_mode = "free";
    end
    if mode ~= "opt" && mode ~= "reach"
        error('RPO2DS_transfer:InvalidMode', "mode must be 'opt' or 'reach'.");
    end
    if params.time_variable_mode ~= "node_times" && params.time_variable_mode ~= "propagation_times"
        error('RPO2DS_transfer:InvalidTimeVariableMode', ...
            "params.time_variable_mode must be 'node_times' or 'propagation_times'.");
    end
    if params.geometry_mode ~= "free_geometry" && params.geometry_mode ~= "prescribed_geometry"
        error('RPO2DS_transfer:InvalidGeometryMode', ...
            "params.geometry_mode must be 'free_geometry' or 'prescribed_geometry'.");
    end
if mode == "reach" && (~isfield(opts, 'dv_budget') || isempty(opts.dv_budget))
    error('RPO2DS_transfer:MissingDvBudget', 'opts.dv_budget is required for reachability mode.');
end
if opts.jacobian_backend ~= "dense" && opts.jacobian_backend ~= "sparse"
    error('RPO2DS_transfer:InvalidJacobianBackend', ...
        "opts.jacobian_backend must be 'dense' or 'sparse'.");
end

% Fixed-arrival timing is parameterized by segment durations anchored at t(N).
if params.operational_mode == "fixed_arrival"
    params.time_variable_mode = "propagation_times";
end

% Build the nominal ECI target used by patched-conic seed generation.
local_progress_log(opts.progress_logs, 'RPO2DS', 'building seed target for target frame %s', char(params.target_frame));
vInfPlusEarthSeed = local_build_seed_vinf_target(vInfPlusEarth, params);
local_progress_log(opts.progress_logs, 'RPO2DS', 'seed target built; |v_inf| %.6f km/s', norm(vInfPlusEarthSeed));

% Reuse the previous NLP state directly when the caller requests continuation.
if local_has_continuation_base(opts)
    local_progress_log(opts.progress_logs, 'RPO2DS', 'continuation path selected');
    result = local_run_continuation_case(mode, vInfPlusEarthSeed, xx0_RPO, params, opts);
    local_progress_log(opts.progress_logs, 'RPO2DS', 'continuation path complete');
    return
end

% Constants (shared)
local_progress_log(opts.progress_logs, 'RPO2DS', 'computing constants');
constants = compute_constants(params);
constants.xx0RPOBase = [xx0_RPO(1); xx0_RPO(2); 0; xx0_RPO(3); xx0_RPO(4); 0]; % Keep the unphased planar RPO seed.
local_progress_log(opts.progress_logs, 'RPO2DS', 'constants ready; TU %.6f days, VU %.6f km/s', constants.TU_days, constants.VU);

% Build the shared phase-policy setup and the database subset used by the solver.
operational_guess = struct();
if params.operational_mode == "standard"
    local_progress_log(opts.progress_logs, 'RPO2DS', 'standard path selected; building phase policy');
    phase_policy = RPO2DS_phase_policy(RPOMoonDatabase, params);
    RPOMoonDatabaseSolve = phase_policy.database;
    local_progress_log(opts.progress_logs, 'RPO2DS', 'phase policy ready; solve database rows %d', height(RPOMoonDatabaseSolve));

    % Find the initial guess
    % Moon-Earth leg parameters
    params_guess = struct('N', params.N, ...
        'hStart', params.hStart, ...
        'sweep_rad', params.sweep_rad, ...
        'nSol', params.nSol, ...
        'feasibleLimit', params.feasibleLimit, ...
        'nMoonOrbitSamples', params.nMoonOrbitSamples, ...
        'time_mode', params.time_mode, ...
        'maxTOF', params.maxTOF, ...
        'maxTransferTOF', params.maxTransferTOF, ...
        'geometry_mode', params.geometry_mode, ...
        'initial_geometry_time', params.initial_geometry_time, ...
        'target_frame', params.target_frame, ...
        'vinf_target_rth_km_s', params.vinf_target_rth_km_s, ...
        'theta_sun0_mode', params.theta_sun0_mode, ...
        'theta_sun0_value', params.theta_sun0_value, ...
        'theta_sun0_seed', params.theta_sun0_seed, ...
        'ws_sun_em_nd', params.ws_sun_em_nd, ...
        'rth_seed_target_iterations', params.rth_seed_target_iterations, ...
        'rth_seed_nominal_lambert_tof_days', params.rth_seed_nominal_lambert_tof_days, ...
        'rth_seed_max_lambert_tof_days', params.rth_seed_max_lambert_tof_days, ...
        'progress_logs', opts.progress_logs);
    local_progress_log(opts.progress_logs, 'RPO2DS', 'calling RPO_Moon_Earth_guess for initial guess');
    [FullSolution, ~, h, guess_diagnostics] = RPO_Moon_Earth_guess(vInfPlusEarthSeed, RPOMoonDatabaseSolve, params_guess);
    local_progress_log(opts.progress_logs, 'RPO2DS', 'initial guess patcher returned; rows %d, h %.6f km', height(FullSolution), h);
    phase_policy = finalize_phase_policy(phase_policy, FullSolution);
    local_progress_log(opts.progress_logs, 'RPO2DS', 'phase policy finalized; seed phase %.6f', phase_policy.phase_seed);
    if ~local_uses_standard_rth_time_aware_seed(params)
        local_progress_log(opts.progress_logs, 'RPO2DS', 'refining standard RTH patch seed');
        [FullSolution, h, guess_diagnostics, vInfPlusEarthSeed, phase_policy] = ...
            local_refine_standard_rth_patch_seed(FullSolution, h, guess_diagnostics, vInfPlusEarthSeed, ...
            RPOMoonDatabaseSolve, xx0_RPO, params, params_guess, constants, opts, phase_policy);
        local_progress_log(opts.progress_logs, 'RPO2DS', 'standard RTH patch seed refinement complete');
    end
else
    local_progress_log(opts.progress_logs, 'RPO2DS', 'operational path selected: %s', char(params.operational_mode));
    phase_template = RPO2DS_phase_policy(RPOMoonDatabase, struct('phase_mode', 'no_phase', 'phase_value', []));
    params.TU_days_operational = constants.TU_days;
    params.progress_logs = opts.progress_logs;
    [FullSolution, ~, h, guess_diagnostics, phase_policy, operational_guess] = ...
        RPO2DS_operational_guess(vInfPlusEarthSeed, RPOMoonDatabase, params, phase_template, constants);
    RPOMoonDatabaseSolve = RPOMoonDatabase;
    local_progress_log(opts.progress_logs, 'RPO2DS', 'operational initial guess returned; rows %d, h %.6f km', height(FullSolution), h);
end
if opts.print_logs
    fprintf('Phase mode: %s\n', char(phase_policy.mode));
    fprintf('Seed phase: %.6f\n', phase_policy.phase_seed);
    if ~isempty(phase_policy.phase_request_wrapped)
        fprintf('Requested phase: %.6f\n', phase_policy.phase_request_wrapped);
    end
end

% Build initial trajectory (shared)
local_progress_log(opts.progress_logs, 'RPO2DS', 'assembling initial trajectory');
traj = assemble_trajectory(FullSolution, xx0_RPO, RPOMoonDatabaseSolve, constants, opts);
traj = local_anchor_operational_trajectory_time(traj, params, phase_policy);
constants.xx0RPO = traj.xx0RPO;
seed_time_info = build_seed_time_info(FullSolution, traj, params, guess_diagnostics, constants);
phase_policy = local_update_operational_phase_seed(phase_policy, params, traj, constants);
local_progress_log(opts.progress_logs, 'RPO2DS', 'initial trajectory ready; seed TOF %.6f days', seed_time_info.seed_tof_exact_days);
if opts.print_logs
    fprintf('Time mode: %s\n', seed_time_info.mode);
    fprintf('Seed exact TOF: %.6f days\n', seed_time_info.seed_tof_exact_days);
    if isfinite(seed_time_info.seed_tof_lower_bound_days)
        fprintf('Seed lower-bound TOF: %.6f days\n', seed_time_info.seed_tof_lower_bound_days);
    end
    fprintf('Requested max TOF: %.6f days\n', seed_time_info.seed_tof_cap_days);
    fprintf('Seed satisfies cap: %d\n', seed_time_info.seed_tof_keep);
    fprintf('Active optimization max TOF: %.6f days\n', seed_time_info.active_max_tof_days);
end

% Node generation and initial guess (shared)
local_progress_log(opts.progress_logs, 'RPO2DS', 'building multiple-shooting nodes');
nodes = build_nodes(traj, params, constants, opts);
[nodes, FullSolution] = local_align_standard_rth_final_node(nodes, FullSolution, params, constants);
local_progress_log(opts.progress_logs, 'RPO2DS', 'nodes ready; final node count %d', nodes.n_nodes);

% Initial guess for regularized impulses (shared)
local_progress_log(opts.progress_logs, 'RPO2DS', 'initializing regularized impulses');
impulses = init_impulses(FullSolution, nodes, constants, params);
local_progress_log(opts.progress_logs, 'RPO2DS', 'regularized impulses ready');

% Pack decision vector and build constraints constants
local_progress_log(opts.progress_logs, 'RPO2DS', 'packing decision vector and NLP data');
[y0, problem_data, meta, lb, ub] = build_setup_data( ...
    mode, vInfPlusEarthSeed, FullSolution, nodes, impulses, constants, params, opts, phase_policy, seed_time_info);
meta = local_append_operational_meta(meta, problem_data, operational_guess);
local_progress_log(opts.progress_logs, 'RPO2DS', 'NLP data packed; decision variables %d, constraints mode %s', numel(y0), problem_data.mode);

% Save the setup package used by tests, scripts, and optional setup-only runs.
setup_data = struct();
setup_data.y0_seed = y0;
setup_data.problem_data = problem_data;
setup_data.lb = lb;
setup_data.ub = ub;

% Replace the internally generated decision-vector seed when the caller provides a warm start.
y0 = local_apply_y0_override(y0, opts);
setup_data.y0 = y0;

% Rebuild the NLP seed so it can be reused by plots and result diagnostics.
local_progress_log(opts.progress_logs, 'RPO2DS', 'rebuilding pre-optimization trajectory from y0');
preopt_results = struct();
preopt_results.setup = local_build_preopt_setup(mode, y0, nodes, impulses, constants, meta, opts);
preopt_results.stage1 = [];
local_progress_log(opts.progress_logs, 'RPO2DS', 'pre-optimization trajectory ready');

% Build the pre-solve plots that are available directly from the initial setup.
plot_flags = opts.graph;
plot_results = struct();
if plot_flags.flyby
    plot_results.flyby = plot_flyby(traj, constants);
end
if plot_flags.initial_guess
    plot_results.initial_guess = plot_initial_guess(traj, constants);
end
if plot_flags.nodes
    plot_results.nodes = plot_nodes(traj, nodes, constants);
end
if plot_flags.initial_guess_fmincon_rot
    preopt_results.setup = local_ensure_full_preopt_post(preopt_results.setup, mode, y0, nodes, impulses, constants, meta, opts);
    plot_results.initial_guess_fmincon_rot = plot_initial_guess_fmincon_rot(preopt_results.setup, constants);
end
if plot_flags.initial_guess_fmincon_inertial
    preopt_results.setup = local_ensure_full_preopt_post(preopt_results.setup, mode, y0, nodes, impulses, constants, meta, opts);
    plot_results.initial_guess_fmincon_inertial = plot_initial_guess_fmincon_inertial(preopt_results.setup, constants);
end

% Stop after the setup stage when the caller only wants the initial architecture.
if opts.setup_only
    result = build_result_struct(mode, FullSolution, vInfPlusEarth, xx0_RPO, params, opts, constants, traj, nodes, impulses, meta, ...
        [], [], preopt_results, plot_results, [], [], h, guess_diagnostics, seed_time_info, [], setup_data);
    return
end

% Print the selected patch and initial-guess diagnostics before the optimization starts.
if opts.print_patch_logs
    local_print_preoptimization_logs(FullSolution, guess_diagnostics, seed_time_info, traj, constants, opts);
end

% Optional gradient check
if opts.gradient_check
    if phase_policy.use_phase_variable
        grad_results = run_gradient_checks_phase(mode, y0, problem_data, opts);
    else
        grad_results = run_gradient_checks(mode, y0, problem_data, opts);
    end
    result = build_result_struct(mode, FullSolution, vInfPlusEarthSeed, xx0_RPO, params, opts, constants, traj, nodes, impulses, meta, ...
        [], [], preopt_results, struct(), [], [], h, guess_diagnostics, seed_time_info, grad_results, setup_data);
    return
end

% Solve (fmincon)
local_progress_log(opts.progress_logs, 'RPO2DS', 'dispatching solver');
[solve_out] = run_solver_dispatch(mode, y0, problem_data, opts, phase_policy, lb, ub, nodes, impulses, constants, meta);
local_progress_log(opts.progress_logs, 'RPO2DS', 'solver returned; exitflag %g', solve_out.exitflag);

% Post-process
local_progress_log(opts.progress_logs, 'RPO2DS', 'post-processing optimized solution');
if phase_policy.use_phase_variable
    post = postprocess_solution_phase(mode, solve_out.y_opt, nodes, impulses, constants, meta, opts);
else
    post = postprocess_solution(mode, solve_out.y_opt, nodes, impulses, constants, meta, opts);
end
if mode == "reach"
    local_progress_log(opts.progress_logs, 'RPO2DS', 'post-processing complete; dv %.6f m/s, tof %.6f days, vInfMag %.6f km/s', post.dv_opt, post.tof_opt, post.vInfMag);
else
    local_progress_log(opts.progress_logs, 'RPO2DS', 'post-processing complete; dv %.6f m/s, tof %.6f days', post.dv_opt, post.tof_opt);
end

% Rebuild the stage-1 feasibility solution for plots and diagnostics.
if isfield(solve_out, 'y_stage1') && ~isempty(solve_out.y_stage1)
    preopt_results.stage1 = build_preopt_post(mode, solve_out.y_stage1, nodes, impulses, constants, meta, opts);
end

% Plots
if isfield(solve_out, 'stage1_plot_results') && ~isempty(solve_out.stage1_plot_results)
    if isfield(solve_out.stage1_plot_results, 'objzero_rot')
        plot_results.objzero_rot = solve_out.stage1_plot_results.objzero_rot;
    end
    if isfield(solve_out.stage1_plot_results, 'objzero_inertial')
        plot_results.objzero_inertial = solve_out.stage1_plot_results.objzero_inertial;
    end
end
if plot_flags.optimized_rot
    plot_results.optimized_rot = plot_optimized_rot(post, constants);
end
if plot_flags.optimized_inertial
    plot_results.optimized_inertial = plot_optimized_inertial(post, constants);
end

% Build result struct
local_progress_log(opts.progress_logs, 'RPO2DS', 'building result struct');
result = build_result_struct(mode, FullSolution, vInfPlusEarth, xx0_RPO, params, opts, constants, traj, nodes, impulses, meta, ...
    solve_out, post, preopt_results, plot_results, solve_out.exitflag, solve_out.output, h, guess_diagnostics, seed_time_info, [], setup_data);
local_progress_log(opts.progress_logs, 'RPO2DS', 'result struct ready');

end

function [params, opts] = apply_defaults(params, opts)
%APPLY_DEFAULTS Normalize the supported parameter forms and fill current solver defaults.
%
% PROTOTYPE:
%   [params, opts] = apply_defaults(params, opts)
%
% DESCRIPTION:
%   Preserve caller-supplied nonempty settings while resolving the existing
%   parameter, graph, backend, geometry, timing, phase, and impulse policies.
%
% INPUT:
%   params  Transfer policy/settings struct; see the public RPO2DS_transfer
%           header.
%   opts  Backend, iteration, tolerance, logging, plotting, and warm-start
%         options.
%
% OUTPUT:
%   params  Normalized transfer policy/settings struct.
%   opts  Normalized backend, graph, and solve options.

if isnumeric(params)
    params_vec = params;
    params = struct();
    if numel(params_vec) >= 1
        params.n_nodes = params_vec(1);
    end
    if numel(params_vec) >= 2
        params.hStart = params_vec(2);
    end
end

% Defaults for params
if ~isfield(params, 'n_nodes') || isempty(params.n_nodes)
    params.n_nodes = 40;
end
if ~isfield(params, 'hStart') || isempty(params.hStart)
    params.hStart = 1000; % km
end
if ~isfield(params, 'N') || isempty(params.N)
    params.N = 101;
end
if ~isfield(params, 'sweep_rad') || isempty(params.sweep_rad)
    params.sweep_rad = 2 * pi;
end
if ~isfield(params, 'nSol') || isempty(params.nSol)
    params.nSol = 1;
end
if ~isfield(params, 'feasibleLimit') || isempty(params.feasibleLimit)
    params.feasibleLimit = 1e-2;
end
if ~isfield(params, 'maxTOF') || isempty(params.maxTOF)
    params.maxTOF = 100; % days
end
if ~isfield(params, 'maxTransferTOF') || isempty(params.maxTransferTOF)
    params.maxTransferTOF = params.maxTOF;
end
if ~isfield(params, 'max_final_earth_distance_km')
    params.max_final_earth_distance_km = [];
end
if ~isfield(params, 'impulse_mode') || isempty(params.impulse_mode)
    params.impulse_mode = 'every_node';
end
if ~isfield(params, 'nMoonOrbitSamples') || isempty(params.nMoonOrbitSamples)
    params.nMoonOrbitSamples = 100;
end
if ~isfield(params, 'operational_mode') || isempty(params.operational_mode)
    params.operational_mode = 'standard';
end
if ~isfield(params, 'time_mode') || isempty(params.time_mode)
    params.time_mode = 'optimization_only';
end
if ~isfield(params, 'maxTOF_relax_factor') || isempty(params.maxTOF_relax_factor)
    params.maxTOF_relax_factor = 1.1;
end
if ~isfield(params, 'time_variable_mode') || isempty(params.time_variable_mode)
    params.time_variable_mode = 'node_times';
end
if ~isfield(params, 'geometry_mode') || isempty(params.geometry_mode)
    params.geometry_mode = 'free_geometry';
end
if ~isfield(params, 'initial_geometry_time') || isempty(params.initial_geometry_time)
    params.initial_geometry_time = 0;
end
if ~isfield(params, 't_node1_fixed')
    params.t_node1_fixed = [];
end
if ~isfield(params, 't_nodeN_fixed')
    params.t_nodeN_fixed = [];
end
if ~isfield(params, 'dense_nodes_multiplier') || isempty(params.dense_nodes_multiplier)
    params.dense_nodes_multiplier = 0.15;
end
if ~isfield(params, 'node_count_mode') || isempty(params.node_count_mode)
    params.node_count_mode = 'fixed';
end
if ~isfield(params, 'nodes_per_day') || isempty(params.nodes_per_day)
    params.nodes_per_day = 0.5;
end
if ~isfield(params, 'n_nodes_min') || isempty(params.n_nodes_min)
    params.n_nodes_min = 10;
end
if ~isfield(params, 'n_nodes_max') || isempty(params.n_nodes_max)
    params.n_nodes_max = 200;
end
if ~isfield(params, 'adaptive_plus_divisions') || isempty(params.adaptive_plus_divisions)
    params.adaptive_plus_divisions = 6;
end
if ~isfield(params, 'adaptive_plus_upstream_nodes') || isempty(params.adaptive_plus_upstream_nodes)
    params.adaptive_plus_upstream_nodes = 2;
end
if ~isfield(params, 'phase_mode') || isempty(params.phase_mode)
    params.phase_mode = 'seed_phase';
end
if ~isfield(params, 'phase_value') || isempty(params.phase_value)
    params.phase_value = [];
end
if ~isfield(params, 'initial_epoch_time') || isempty(params.initial_epoch_time)
    params.initial_epoch_time = 0;
end
if ~isfield(params, 'initial_rpo_phase') || isempty(params.initial_rpo_phase)
    params.initial_rpo_phase = 0;
end
if ~isfield(params, 'final_epoch_time')
    params.final_epoch_time = [];
end
if ~isfield(params, 'operational_seed_pool_size') || isempty(params.operational_seed_pool_size)
    params.operational_seed_pool_size = max(params.nSol, 12);
end
if ~isfield(params, 'fixed_arrival_phase_tolerance') || isempty(params.fixed_arrival_phase_tolerance)
    params.fixed_arrival_phase_tolerance = [];
end
if strcmpi(string(params.operational_mode), "fixed_arrival") && isempty(params.final_epoch_time)
    params.final_epoch_time = params.initial_epoch_time + params.maxTOF * 86400 / (1 / 2.665e-6);
end
if ~isfield(params, 'target_frame') || isempty(params.target_frame)
    params.target_frame = 'eci';
end
if ~isfield(params, 'vinf_target_rth_km_s') || isempty(params.vinf_target_rth_km_s)
    params.vinf_target_rth_km_s = [0; 3; 0];
end
if ~isfield(params, 'theta_sun0_mode') || isempty(params.theta_sun0_mode)
    params.theta_sun0_mode = 'specified';
end
if ~isfield(params, 'theta_sun0_value') || isempty(params.theta_sun0_value)
    params.theta_sun0_value = 0;
end
if ~isfield(params, 'theta_sun0_seed') || isempty(params.theta_sun0_seed)
    params.theta_sun0_seed = params.theta_sun0_value;
end
if ~isfield(params, 'ws_sun_em_nd') || isempty(params.ws_sun_em_nd)
    params.ws_sun_em_nd = -9.25195985e-1;
end
if ~isfield(params, 'rth_seed_target_iterations') || isempty(params.rth_seed_target_iterations)
    params.rth_seed_target_iterations = 2;
end
if ~isfield(params, 'rth_seed_nominal_lambert_tof_days') || isempty(params.rth_seed_nominal_lambert_tof_days)
    params.rth_seed_nominal_lambert_tof_days = 8;
end
if ~isfield(params, 'rth_seed_max_lambert_tof_days') || isempty(params.rth_seed_max_lambert_tof_days)
    params.rth_seed_max_lambert_tof_days = params.maxTransferTOF;
end
if ~isfield(params, 'final_time_mode') || isempty(params.final_time_mode)
    params.final_time_mode = 'free';
    if strcmpi(string(params.operational_mode), "fixed_arrival")
        params.final_time_mode = 'specified';
    end
end
if ~isfield(params, 'final_time_seed') || isempty(params.final_time_seed)
    if ~isempty(params.final_epoch_time)
        params.final_time_seed = params.final_epoch_time;
    else
        params.final_time_seed = params.initial_epoch_time + params.maxTOF * 86400 / (1 / 2.665e-6);
    end
end

% Defaults for opts
if ~isfield(opts, 'gradient_check') || isempty(opts.gradient_check)
    opts.gradient_check = 0;
end
if ~isfield(opts, 'jacobian_backend') || isempty(opts.jacobian_backend)
    opts.jacobian_backend = "dense";
else
    opts.jacobian_backend = lower(string(opts.jacobian_backend));
end
if ~isfield(opts, 'solver_backend') || isempty(opts.solver_backend)
    opts.solver_backend = "fmincon";
else
    opts.solver_backend = lower(string(opts.solver_backend));
end
if ~isfield(opts, 'casadi_path') || isempty(opts.casadi_path)
    % CasADi is resolved by rpo2ds.require_solver from the MATLAB path or
    % the repository's ignored third_party/casadi directory.  Do not embed a
    % machine-specific installation path in the public backend.
    opts.casadi_path = '';
end
if ~isfield(opts, 'solver_display') || isempty(opts.solver_display)
    opts.solver_display = 'iter';
end
if ~isfield(opts, 'stage1_algorithm') || isempty(opts.stage1_algorithm)
    opts.stage1_algorithm = 'sqp';
end
if ~isfield(opts, 'stage2_algorithm') || isempty(opts.stage2_algorithm)
    opts.stage2_algorithm = 'sqp';
end
if ~isfield(opts, 'max_iter_stage1') || isempty(opts.max_iter_stage1)
    opts.max_iter_stage1 = 2e3;
end
if ~isfield(opts, 'max_iter_stage2') || isempty(opts.max_iter_stage2)
    opts.max_iter_stage2 = 2e3;
end
if ~isfield(opts, 'constraint_tolerance') || isempty(opts.constraint_tolerance)
    opts.constraint_tolerance = 1e-6;
end
if ~isfield(opts, 'optimality_tolerance') || isempty(opts.optimality_tolerance)
    opts.optimality_tolerance = 1e-4;
end
if ~isfield(opts, 'print_logs') || isempty(opts.print_logs)
    opts.print_logs = true;
end
if ~isfield(opts, 'print_patch_logs') || isempty(opts.print_patch_logs)
    opts.print_patch_logs = true;
end
if ~isfield(opts, 'progress_logs') || isempty(opts.progress_logs)
    opts.progress_logs = false;
end
if ~isfield(opts, 'patch_log_candidate_rows') || isempty(opts.patch_log_candidate_rows)
    opts.patch_log_candidate_rows = inf;
end
if ~isfield(opts, 'setup_only') || isempty(opts.setup_only)
    opts.setup_only = false;
end

% Normalize graph flags
if ~isfield(opts, 'graph') || isempty(opts.graph)
    opts.graph = struct();
end
opts.graph = normalize_graph_flags(opts.graph);

end

function [mode, params, opts] = local_apply_casadi_ipopt_profile_v1(mode_in, params, opts)
%LOCAL_APPLY_CASADI_IPOPT_PROFILE_V1 Normalize settings onto the supported CasADi/IPOPT profile.
%
% PROTOTYPE:
%   [mode, params, opts] = local_apply_casadi_ipopt_profile_v1(mode_in,
%   params, opts)
%
% DESCRIPTION:
%   Validate supported mode/phase/time combinations before building the symbolic
%   NLP. Preserve the selected target frame and operational timing controls.
%
% INPUT:
%   mode_in  Requested opt/reach mode before applying the backend profile.
%   params  Transfer policy/settings struct; see the public RPO2DS_transfer
%           header.
%   opts  Backend, iteration, tolerance, logging, plotting, and warm-start
%         options.
%
% OUTPUT:
%   mode  Normalized opt/reach mode after backend-profile selection.
%   params  Normalized transfer policy/settings struct.
%   opts  Normalized backend, graph, and solve options.

requested_mode = lower(string(mode_in));
requested_operational_mode = lower(string(params.operational_mode));
requested_phase_mode = lower(string(params.phase_mode));
requested_geometry_mode = lower(string(params.geometry_mode));

if opts.print_logs
    fprintf('IPOPT mode selected. Normalizing settings for CasADi profile v1:\n');
    fprintf('  mode kept as: %s\n', char(requested_mode));
    fprintf('  params.time_variable_mode set to: propagation_times\n');
    fprintf('  params.time_mode set to: max_time\n');
    fprintf('  opts.jacobian_backend ignored: CasADi automatic differentiation is used\n');
    fprintf('  opts.stage1_algorithm/stage2_algorithm ignored: IPOPT is used\n');
end

% Keep the requested solve family now that the CasADi builder supports reach mode.
mode = requested_mode;
params.time_variable_mode = 'propagation_times';
params.time_mode = 'max_time';
opts.jacobian_backend = "dense";

% Preserve the standard departure-driven profiles now covered by the CasADi builder.
if requested_operational_mode == "standard"
    if opts.print_logs
        fprintf('  params.operational_mode kept as: standard\n');
        fprintf('  params.phase_mode kept as: %s\n', char(requested_phase_mode));
        fprintf('  params.geometry_mode kept as: %s\n', char(requested_geometry_mode));
    end
    params.operational_mode = 'standard';
    params.phase_mode = char(requested_phase_mode);
    params.geometry_mode = char(requested_geometry_mode);
    opts.casadi_ipopt_profile = "standard";
    if opts.print_logs && requested_phase_mode == "other_phase"
        fprintf('  phase_value is kept fixed by NLP bounds\n');
    end
    if opts.print_logs && requested_geometry_mode == "prescribed_geometry"
        fprintf('  initial_geometry_time is kept fixed by NLP bounds\n');
    end
    return
end

% Preserve the operational fixed-arrival timing profile when it is requested.
if requested_operational_mode == "fixed_arrival"
    if opts.print_logs
        fprintf('  params.operational_mode kept as: fixed_arrival\n');
        fprintf('  params.geometry_mode ignored: fixed_arrival controls the timing geometry\n');
        fprintf('  params.phase_mode ignored: departure phase is tied to arrival epoch\n');
    end
    params.operational_mode = 'fixed_arrival';
    params.geometry_mode = 'free_geometry';
    params.phase_mode = 'no_phase';
    params.phase_value = [];
    opts.casadi_ipopt_profile = "fixed_arrival";
    return
end

% Otherwise keep the original initial-state IPOPT profile.
if opts.print_logs
    fprintf('  params.operational_mode set to: initial_state\n');
    fprintf('  params.geometry_mode ignored: initial_state controls the timing geometry\n');
    fprintf('  params.phase_mode ignored: departure phase is derived from initial_epoch_time and initial_rpo_phase\n');
end
params.operational_mode = 'initial_state';
params.geometry_mode = 'free_geometry';
params.phase_mode = 'no_phase';
params.phase_value = [];
opts.casadi_ipopt_profile = "initial_state";
end

function has_base = local_has_continuation_base(opts)
%LOCAL_HAS_CONTINUATION_BASE Check whether the caller provided one previous solver result for continuation.

has_base = isfield(opts, 'continuation_base_result') && ~isempty(opts.continuation_base_result);
end

function result = local_run_continuation_case(mode, vInfPlusEarth, xx0_RPO, params, opts)
%LOCAL_RUN_CONTINUATION_CASE Reuse a solved NLP state for a new target or timing constraint.
%
% PROTOTYPE:
%   result = local_run_continuation_case(mode, vInfPlusEarth, xx0_RPO,
%   params, opts)
%
% DESCRIPTION:
%   Rebuild the continuation setup from the previous result, update active
%   target/epoch fields, then use the shared solve and postprocessing path.
%
% INPUT:
%   mode  opt minimizes regularized delta-v; reach maximizes v-infinity.
%   vInfPlusEarth  Earth-relative ECI target velocity/direction [km/s].
%   xx0_RPO  Planar rotating seed [x; y; vx; vy] in DU/VU.
%   params  Transfer policy/settings struct; see the public RPO2DS_transfer
%           header.
%   opts  Backend, iteration, tolerance, logging, plotting, and warm-start
%         options.
%
% OUTPUT:
%   result  Full transfer result with inputs, setup/solve provenance,
%           trajectories, and metrics.

local_progress_log(opts.progress_logs, 'RPO2DS', 'continuation setup: validating base result');
base_result = opts.continuation_base_result;
if ~isstruct(base_result) || ~isfield(base_result, 'setup') || ~isfield(base_result, 'meta')
    error('RPO2DS_transfer:InvalidContinuationBase', ...
        'opts.continuation_base_result must be one previous RPO2DS_transfer result struct.');
end

% Reuse the previous shared structures so continuation feeds directly into the next fmincon cycle.
constants = base_result.constants;
traj = base_result.initial_guess;
nodes = base_result.nodes;
impulses = base_result.impulses;
FullSolution = base_result.FullSolution;
h = base_result.h;
guess_diagnostics = base_result.initial_guess_search;
phase_policy = struct('use_phase_variable', logical(base_result.meta.use_phase_variable));

% Start from the previous optimized solution and only update the fixed initial-time constraint.
local_progress_log(opts.progress_logs, 'RPO2DS', 'continuation setup: rebuilding y0/problem/bounds');
[y0, problem_data, meta, lb, ub, setup_data] = local_build_continuation_setup(base_result, vInfPlusEarth, params, opts);
local_progress_log(opts.progress_logs, 'RPO2DS', 'continuation setup ready; decision variables %d', numel(y0));

% Rebuild the NLP seed so it can be reused by plots and result diagnostics.
preopt_results = struct();
local_progress_log(opts.progress_logs, 'RPO2DS', 'continuation setup: rebuilding preopt trajectory');
preopt_results.setup = local_build_preopt_setup(mode, y0, nodes, impulses, constants, meta, opts);
preopt_results.stage1 = [];
local_progress_log(opts.progress_logs, 'RPO2DS', 'continuation preopt trajectory ready');

% Build the pre-solve plots that are available directly from the continuation setup.
plot_flags = opts.graph;
plot_results = struct();
if plot_flags.flyby
    plot_results.flyby = plot_flyby(traj, constants);
end
if plot_flags.initial_guess
    plot_results.initial_guess = plot_initial_guess(traj, constants);
end
if plot_flags.nodes
    plot_results.nodes = plot_nodes(traj, nodes, constants);
end
if plot_flags.initial_guess_fmincon_rot
    preopt_results.setup = local_ensure_full_preopt_post(preopt_results.setup, mode, y0, nodes, impulses, constants, meta, opts);
    plot_results.initial_guess_fmincon_rot = plot_initial_guess_fmincon_rot(preopt_results.setup, constants);
end
if plot_flags.initial_guess_fmincon_inertial
    preopt_results.setup = local_ensure_full_preopt_post(preopt_results.setup, mode, y0, nodes, impulses, constants, meta, opts);
    plot_results.initial_guess_fmincon_inertial = plot_initial_guess_fmincon_inertial(preopt_results.setup, constants);
end

% Stop after rebuilding the continuation NLP setup when the caller only wants the architecture.
seed_time_info = local_build_seed_time_info_from_result(base_result);
if opts.setup_only
    result = build_result_struct(mode, FullSolution, vInfPlusEarth, xx0_RPO, params, opts, constants, traj, nodes, impulses, meta, ...
        [], [], preopt_results, plot_results, [], [], h, guess_diagnostics, seed_time_info, [], setup_data);
    return
end

% Print the selected patch and the reused continuation setup when compact logs are enabled.
if opts.print_patch_logs
    local_print_preoptimization_logs(FullSolution, guess_diagnostics, seed_time_info, traj, constants, opts);
end

% Run the shared solver path directly on the reused NLP state.
local_progress_log(opts.progress_logs, 'RPO2DS', 'continuation dispatching solver');
[solve_out] = run_solver_dispatch(mode, y0, problem_data, opts, phase_policy, lb, ub, nodes, impulses, constants, meta);
local_progress_log(opts.progress_logs, 'RPO2DS', 'continuation solver returned; exitflag %g', solve_out.exitflag);

% Post-process the continuation solution with the reused shared geometry structs.
local_progress_log(opts.progress_logs, 'RPO2DS', 'continuation post-processing optimized solution');
if phase_policy.use_phase_variable
    post = postprocess_solution_phase(mode, solve_out.y_opt, nodes, impulses, constants, meta, opts);
else
    post = postprocess_solution(mode, solve_out.y_opt, nodes, impulses, constants, meta, opts);
end
if mode == "reach"
    local_progress_log(opts.progress_logs, 'RPO2DS', 'continuation post-processing complete; dv %.6f m/s, tof %.6f days, vInfMag %.6f km/s', post.dv_opt, post.tof_opt, post.vInfMag);
else
    local_progress_log(opts.progress_logs, 'RPO2DS', 'continuation post-processing complete; dv %.6f m/s, tof %.6f days', post.dv_opt, post.tof_opt);
end

% Rebuild the stage-1 feasibility solution for plots and diagnostics.
if isfield(solve_out, 'y_stage1') && ~isempty(solve_out.y_stage1)
    preopt_results.stage1 = build_preopt_post(mode, solve_out.y_stage1, nodes, impulses, constants, meta, opts);
end

% Rebuild the shared plots for the continuation solve.
if isfield(solve_out, 'stage1_plot_results') && ~isempty(solve_out.stage1_plot_results)
    if isfield(solve_out.stage1_plot_results, 'objzero_rot')
        plot_results.objzero_rot = solve_out.stage1_plot_results.objzero_rot;
    end
    if isfield(solve_out.stage1_plot_results, 'objzero_inertial')
        plot_results.objzero_inertial = solve_out.stage1_plot_results.objzero_inertial;
    end
end
if plot_flags.optimized_rot
    plot_results.optimized_rot = plot_optimized_rot(post, constants);
end
if plot_flags.optimized_inertial
    plot_results.optimized_inertial = plot_optimized_inertial(post, constants);
end

% Build the shared result struct using the updated continuation metadata and the reused seed diagnostics.
local_progress_log(opts.progress_logs, 'RPO2DS', 'continuation building result struct');
result = build_result_struct(mode, FullSolution, vInfPlusEarth, xx0_RPO, params, opts, constants, traj, nodes, impulses, meta, ...
    solve_out, post, preopt_results, plot_results, solve_out.exitflag, solve_out.output, h, guess_diagnostics, seed_time_info, [], setup_data);
local_progress_log(opts.progress_logs, 'RPO2DS', 'continuation result struct ready');
end

function graph = normalize_graph_flags(graph_in)
%NORMALIZE_GRAPH_FLAGS Resolve the complete supported set of plotting switches.
%
% PROTOTYPE:
%   graph = normalize_graph_flags(graph_in)
%
% DESCRIPTION:
%   Accept the existing flag or struct interface and fill missing plot toggles.
%   These options control diagnostic figures without changing the NLP layout.
%
% INPUT:
%   graph_in  Plotting flag/struct accepted by the public solver interface.
%
% OUTPUT:
%   graph  Complete scalar struct of supported plotting flags.

base = struct('flyby', false, ...
    'initial_guess', false, ...
    'nodes', false, ...
    'initial_guess_fmincon_rot', false, ...
    'initial_guess_fmincon_inertial', false, ...
    'objzero_rot', false, ...
    'objzero_inertial', false, ...
    'optimized_rot', false, ...
    'optimized_inertial', false);

if islogical(graph_in) || isnumeric(graph_in)
    onoff = logical(graph_in);
    graph = base;
    f = fieldnames(graph);
    for k = 1:numel(f)
        graph.(f{k}) = onoff;
    end
    return
end

if ~isstruct(graph_in)
    graph = base;
    return
end

graph = base;
fields = fieldnames(graph_in);
for k = 1:numel(fields)
    graph.(fields{k}) = logical(graph_in.(fields{k}));
end
end

function constants = compute_constants(params)
%COMPUTE_CONSTANTS Collect frozen physical constants and DU/TU/VU scales.
%
% PROTOTYPE:
%   constants = compute_constants(params)
%
% DESCRIPTION:
%   Preserve the transfer method's historical Earth-Moon constants and minimum
%   altitude conversions; these differ slightly from run01 normalization.
%
% INPUT:
%   params  Transfer policy/settings struct; see the public RPO2DS_transfer
%           header.
%
% OUTPUT:
%   constants  Physical constants and the historical DU/TU/VU normalization.

constants.mu_EM = 0.01215;              % Earth-Moon system mass parameter [nondimensional]
constants.muE = 398600;                 % Earth gravitational parameter [km^3/s^2]
constants.DU = 384400;                  % Distance between Earth and Moon [km]
constants.om_EM = 2.665e-6;             % Angular velocity of Earth and Moon around EM barycenter [rad/s]
constants.TU = 1 / constants.om_EM;     % Time unit in Earth-Moon system [s]
constants.TU_days = constants.TU / 86400; % Time unit in days
constants.VU = constants.DU / constants.TU; % Velocity unit [km/s]
constants.Re = 6378;                    % Earth radius [km]
constants.Rm = 1738;                    % Moon radius [km]
constants.rSOI = (constants.mu_EM / (1 - constants.mu_EM))^(2/5); % Moon SOI radius [DU]
constants.hminEarth = params.hStart;    % Min altitude from Earth [km]
constants.hminMoon = 100;               % Min altitude from Moon [km]
constants.rminEarth = constants.Re + constants.hminEarth;
constants.rminMoon = constants.Rm + constants.hminMoon;
constants.rminEarth_DU = constants.rminEarth / constants.DU;
constants.rminMoon_DU = constants.rminMoon / constants.DU;
constants.maxTOF = params.maxTOF * 86400 / constants.TU; % [TU]
constants.time_variable_mode = char(lower(string(params.time_variable_mode)));
end

function vInfSeed = local_build_seed_vinf_target(vInfPlusEarth, params)
%LOCAL_BUILD_SEED_VINF_TARGET Build the ECI target used by the patched-conic initial guess.

% Keep the current inertial target untouched in the default mode.
if lower(string(params.target_frame)) == "eci"
    vInfSeed = vInfPlusEarth(:);
    return
end

% Select the nominal Sun phase used before the NLP chooses its variables.
if lower(string(params.theta_sun0_mode)) == "free"
    theta_sun0_seed = params.theta_sun0_seed;
else
    theta_sun0_seed = params.theta_sun0_value;
end

% Convert the nominal Earth-Sun RTH target into the ECI seed target.
theta_sun_final = theta_sun0_seed + params.ws_sun_em_nd * params.final_time_seed;
vInfSeed = earth_sun_rth_to_eci(params.vinf_target_rth_km_s, theta_sun_final, params.final_time_seed);
end

function uses_time_aware = local_uses_standard_rth_time_aware_seed(params)
%LOCAL_USES_STANDARD_RTH_TIME_AWARE_SEED Check whether standard RTH seeds were already built at candidate epochs.

uses_time_aware = strcmpi(string(params.operational_mode), "standard") && ...
    strcmpi(string(params.target_frame), "earth_sun_rth");
end

function has_specified = local_has_specified_final_time(params)
%LOCAL_HAS_SPECIFIED_FINAL_TIME Check whether the user requested one fixed Earth-perigee epoch.

has_specified = isfield(params, 'final_time_mode') && ...
    strcmpi(string(params.final_time_mode), "specified") && ...
    isfield(params, 'final_time_seed') && ~isempty(params.final_time_seed);
end

function [FullSolution, h, guess_diagnostics, vInfSeed, phase_policy] = local_refine_standard_rth_patch_seed( ...
    FullSolution, h, guess_diagnostics, vInfSeed, RPOMoonDatabaseSolve, xx0_RPO, params, params_guess, constants, opts, phase_policy)
%LOCAL_REFINE_STANDARD_RTH_PATCH_SEED Rebuild the standard patch after the free final epoch is known.

% Only standard Earth-Sun RTH targeting needs this seed-consistency pass.
if ~local_needs_standard_rth_seed_refinement(params)
    return
end

% Use two passes to reduce the mismatch before the final-node repair.
for iter = 1:2
    traj_seed = assemble_trajectory(FullSolution, xx0_RPO, RPOMoonDatabaseSolve, constants, opts);
    traj_seed = local_anchor_operational_trajectory_time(traj_seed, params, phase_policy);
    t_final_seed = local_standard_rth_final_epoch_from_traj(traj_seed, params);
    vInfUpdated = local_build_seed_vinf_target_at_epoch(params, t_final_seed);
    if norm(vInfUpdated - vInfSeed) <= 1e-8
        return
    end
    vInfSeed = vInfUpdated;
    [FullSolution, ~, h, guess_diagnostics] = RPO_Moon_Earth_guess(vInfSeed, RPOMoonDatabaseSolve, params_guess);
    phase_policy = finalize_phase_policy(phase_policy, FullSolution);
end
end

function needs_refinement = local_needs_standard_rth_seed_refinement(params)
%LOCAL_NEEDS_STANDARD_RTH_SEED_REFINEMENT Check whether standard mode has a free final epoch and an RTH target.

needs_refinement = strcmpi(string(params.operational_mode), "standard") && ...
    strcmpi(string(params.target_frame), "earth_sun_rth") && ...
    ~local_has_specified_final_time(params);
end

function vInfSeed = local_build_seed_vinf_target_at_epoch(params, t_final_seed)
%LOCAL_BUILD_SEED_VINF_TARGET_AT_EPOCH Convert the RTH target into ECI at one provisional final epoch.

% Select the Sun phase used by the current seed pass.
if lower(string(params.theta_sun0_mode)) == "free"
    theta_sun0_seed = params.theta_sun0_seed;
else
    theta_sun0_seed = params.theta_sun0_value;
end

% Rotate the RTH target with the provisional final geometry.
theta_sun_final = theta_sun0_seed + params.ws_sun_em_nd * t_final_seed;
vInfSeed = earth_sun_rth_to_eci(params.vinf_target_rth_km_s, theta_sun_final, t_final_seed);
end

function t_final_seed = local_standard_rth_final_epoch_from_traj(traj_seed, params)
%LOCAL_STANDARD_RTH_FINAL_EPOCH_FROM_TRAJ Read the final epoch as it will be represented by the standard NLP layout.

% Prescribed geometry shifts the whole decision timeline to the requested t1.
if strcmpi(string(params.geometry_mode), "prescribed_geometry")
    transfer_time = traj_seed.tTot(traj_seed.earth_peri_idx) - traj_seed.tTot(1);
    t_final_seed = params.initial_geometry_time + transfer_time;
    return
end

% Free-geometry standard mode keeps the raw patched-conic absolute timeline.
t_final_seed = traj_seed.tTot(traj_seed.earth_peri_idx);
end

function traj = assemble_trajectory(FullSolution, xx0_RPO, RPOMoonDatabase, constants, opts)
%ASSEMBLE_TRAJECTORY Reconstruct and join the patched RPO-Moon-Earth trajectory.
%
% PROTOTYPE:
%   traj = assemble_trajectory(FullSolution, xx0_RPO, RPOMoonDatabase,
%   constants, opts)
%
% DESCRIPTION:
%   Build the lunar hyperbolic periapsis state, propagate to the SOI boundaries,
%   and align the RPO and two-body Earth legs in a common rotating timeline.
%   Retain the forced Moon/Earth periapsis indices used by node generation.
%
% INPUT:
%   FullSolution  Selected patch table with physical flyby states, phases,
%                 and times.
%   xx0_RPO  Planar rotating seed [x; y; vx; vy] in DU/VU.
%   RPOMoonDatabase  RPO-Moon candidate table; schema is described in
%                    docs/DATA.md.
%   constants  Shared physical/scaling struct: DU [km], TU [s], VU [km/s],
%              mu_EM.
%   opts  Backend, iteration, tolerance, logging, plotting, and warm-start
%         options.
%
% OUTPUT:
%   traj  Joined RPO, flyby, and Earth-leg arcs in rotating DU/VU with TU
%         times.

mu_EM = constants.mu_EM;
muE = constants.muE;
DU = constants.DU;
VU = constants.VU;
rSOI = constants.rSOI;

% Flyby insertion
v_in = FullSolution.vInfMinusMoon;
v_out = FullSolution.vInfPlusMoon;
posMoon = FullSolution.RI;
[Rp, Vp] = get_perilune_state(v_in, v_out);
RpECI = RTH2ECI_moon(posMoon, Rp, [], 1);
VpECI = RTH2ECI_moon(posMoon, Vp, [], 0);

% Propagate forward and backward the perilune state until it crosses the SOI
x0 = [RpECI(:)/DU; VpECI(:)/VU];
time = FullSolution.trueAnomalyMoonFlyby;
xx0_peri = in2rot(x0', time, mu_EM);

tspanF = [time, time+1];
tspanB = [time, time-1];
opts_ode = odeset('RelTol', 1e-10, 'AbsTol', 1e-12, ...
    'Events', @(t,x) soi_event_CR3BP(t,x,rSOI,mu_EM));

[tF, xxF] = ode113(@(t,x) rhs_CR3BP_3D(t,x,mu_EM), tspanF, xx0_peri, opts_ode);
[tB, xxB] = ode113(@(t,x) rhs_CR3BP_3D(t,x,mu_EM), tspanB, xx0_peri, opts_ode);

xxFB = [xxB(end:-1:2,:); xxF];
ttFB = [tB(end:-1:2,:); tF];

% Total flyby time
traj.tIncFlyby = tB(1) - tB(end);
traj.tTotFlyby = ttFB(end) - ttFB(1);

% CR3BP propagation (RPO-Moon)
xx0RPO = [xx0_RPO(1), xx0_RPO(2), 0, xx0_RPO(3), xx0_RPO(4), 0]';
if ismember("orbitPhase", RPOMoonDatabase.Properties.VariableNames)
    phase = FullSolution.phase;
    xx0RPO = propagate(xx0RPO, 0, phase, mu_EM, '3D', 'CR3BP', 'false');
end

% Adjust time of departure to match the Moon-side handoff epoch.
angle = ttFB(1);
x0 = xx0RPO;
x0(4:6) = x0(4:6) + FullSolution.dv0CR3BP';
if ismember('T_RPOMoonOriginal', FullSolution.Properties.VariableNames)
    angleDifference = angle - FullSolution.T_RPOMoon;
else
    angleDifference = angdiff(FullSolution.trueAnomalyMoonRPO, angle);
end
traj.t0 = angleDifference;
traj.tf = FullSolution.T_RPOMoon + traj.t0;
[~,~,~, xx, tt]  = propagate(x0, traj.t0, traj.tf, mu_EM, '3D', 'CR3BP', 'false');

% True time flyby
traj.ttFB_from_zero = ttFB - ttFB(1);
traj.ttFB_true = traj.ttFB_from_zero + traj.tf;

% TBP propagation (Moon-Earth leg)
rhs = @(t,y) ode_2bp(t,y,muE);
xx0 = [FullSolution.RI'; FullSolution.VI'];
t0_lamb = (tt(end) + traj.tIncFlyby) * constants.TU;
tf_lamb = FullSolution.T_lambert + t0_lamb;
[t, y] = ode113(rhs, [t0_lamb tf_lamb], xx0, opts_ode);

% Keep only after SOI
traj.tSOI = (tt(end) + traj.tTotFlyby);
t = t/constants.TU;
indices = find(t > traj.tSOI);
y = y(indices, :);
t = t(indices);

y_adim(:,1:3) = y(:,1:3) / DU;
y_adim(:,4:6) = y(:,4:6) / VU;
xxRot = in2rot(y_adim, t, mu_EM);

% Put together times and vectors
traj.tTot = [tt; traj.ttFB_true(2:end); t];
traj.xxTot = [xx; xxFB(2:end,:); xxRot];

% Add DVp to last state (for fmincon)
vp = y(end,4:6);
rp = y(end,1:3);
vpaft = vp + FullSolution.deltavp;
xxpaft_rot = in2rot([rp/DU, vpaft/VU], t(end), mu_EM);
traj.xxTot(end,:) = xxpaft_rot;

% Pass to the next phase the index of the RPO-Moon phase
traj.RPO_idx = numel(tt);
traj.moon_peri_idx = traj.RPO_idx + numel(tB) - 2; % -2 for lost indexes in merging
traj.earth_peri_idx = numel(traj.tTot); % The Earth perigee is the last patched state in the assembled transfer

% Check for duplicate values in tTot
[unique_tTot, idx] = unique(traj.tTot);
if opts.print_logs
    if numel(traj.tTot) ~= numel(unique_tTot)
        duplicate_indices = setdiff(1:numel(traj.tTot), unique(idx));
        fprintf('tTot contains duplicate values at indices %s\n', mat2str(duplicate_indices));
    else
        fprintf('tTot does not contain duplicate values.\n');
    end
end

% Store in trajectory
traj.xxF = xxF;
traj.xxB = xxB;
traj.xxFB = xxFB;
traj.ttFB = ttFB;
traj.xx = xx;
traj.tt = tt;
traj.xxRot = xxRot;
traj.t = t;
traj.y = y;
traj.tB = tB;
traj.tF = tF;
traj.xx0RPO = xx0RPO;
if ismember('earthEpochNd', FullSolution.Properties.VariableNames)
    traj.seed_earth_epoch_nd = FullSolution.earthEpochNd(1);
end
if ismember('selectedGeometryEpoch', FullSolution.Properties.VariableNames)
    traj.seed_departure_epoch_nd = FullSolution.selectedGeometryEpoch(1);
end

% Optional plots for flyby and initial guess handled elsewhere
if opts.graph.flyby || opts.graph.initial_guess
    % plots will be done later
end

end

function nodes = build_nodes(traj, params, constants, opts)
%BUILD_NODES Sample shooting nodes while retaining forced periapsis nodes.
%
% PROTOTYPE:
%   nodes = build_nodes(traj, params, constants, opts)
%
% DESCRIPTION:
%   Apply fixed, adaptive, or adaptive-plus timing, preserve the lunar flyby
%   geometry, and interpolate the seed trajectory at the selected epochs.
%
% INPUT:
%   traj  Joined patched trajectory with rotating DU/VU states and TU
%         timelines.
%   params  Transfer policy/settings struct; see the public RPO2DS_transfer
%           header.
%   constants  Shared physical/scaling struct: DU [km], TU [s], VU [km/s],
%              mu_EM.
%   opts  Backend, iteration, tolerance, logging, plotting, and warm-start
%         options.
%
% OUTPUT:
%   nodes  Shooting states/times plus periapsis indices and sampling
%          metadata.

[n_nodes, node_count_info] = local_resolve_node_count(traj, params, constants);

% Build the node times with the selected counting policy.
if strcmp(node_count_info.mode, 'adaptive_plus')
    [t_nodes_IG, moon_peri_idx_nodes, earth_peri_idx_nodes, density_sep_idx, node_count_info] = ...
        local_build_adaptive_plus_node_times(traj, params, constants, node_count_info);
else
    [t_nodes_IG, moon_peri_idx_nodes, earth_peri_idx_nodes, density_sep_idx] = ...
        local_build_standard_node_times(traj, params, n_nodes);
end

% Interpolate the transfer states at the selected node times.
X_nodes_IG = interp1(traj.tTot, traj.xxTot, t_nodes_IG, 'pchip');

% Replace the exact Moon-peri state after interpolation.
X_nodes_IG(moon_peri_idx_nodes, :) = traj.xxTot(traj.moon_peri_idx, :);

% Replace the exact Earth-perigee state after interpolation.
X_nodes_IG(earth_peri_idx_nodes, :) = traj.xxTot(traj.earth_peri_idx, :);

% Update the final node count after the selected construction path.
n_nodes = numel(t_nodes_IG);

nodes = struct();
nodes.n_nodes = n_nodes;
nodes.t_nodes_IG = t_nodes_IG;
nodes.X_nodes_IG = X_nodes_IG;
nodes.moon_peri_idx_nodes = moon_peri_idx_nodes;
nodes.earth_peri_idx_nodes = earth_peri_idx_nodes;
nodes.density_sep_idx = density_sep_idx;
nodes.node_count_mode = node_count_info.mode;
nodes.transfer_time_days = node_count_info.transfer_time_days;
nodes.nodes_per_day = node_count_info.nodes_per_day;
nodes.n_nodes_base = node_count_info.n_nodes_base;
nodes.n_nodes_before_peri = numel(t_nodes_IG) - 1;
nodes.n_nodes_after_peri = n_nodes;
nodes.n_nodes_min = node_count_info.n_nodes_min;
nodes.n_nodes_max = node_count_info.n_nodes_max;
nodes.adaptive_plus_divisions = node_count_info.adaptive_plus_divisions;
nodes.adaptive_plus_upstream_nodes_requested = node_count_info.adaptive_plus_upstream_nodes_requested;
nodes.adaptive_plus_upstream_nodes_used = node_count_info.adaptive_plus_upstream_nodes_used;
nodes.adaptive_plus_dt = node_count_info.adaptive_plus_dt;
nodes.adaptive_plus_dt_days = node_count_info.adaptive_plus_dt_days;
nodes.adaptive_plus_prefix_node_count = node_count_info.adaptive_plus_prefix_node_count;
nodes.adaptive_plus_even_start_idx_nodes = node_count_info.adaptive_plus_even_start_idx_nodes;
nodes.adaptive_plus_between_internal_nodes = node_count_info.adaptive_plus_between_internal_nodes;

if opts.graph.nodes
    % plotting handled elsewhere
end
end

function [nodes, FullSolution] = local_align_standard_rth_final_node(nodes, FullSolution, params, constants)
%LOCAL_ALIGN_STANDARD_RTH_FINAL_NODE Rebuild the final Earth-perigee node for the active free-final-time RTH target.

% Keep the legacy node seed unless the RTH target has a free final epoch.
if ~local_needs_free_final_rth_alignment(params)
    return
end

% Read the final epoch as it will be decoded by the NLP decision vector.
t_final = local_free_final_rth_epoch_from_nodes(nodes, params);
Xn_old_in = rot2in(nodes.X_nodes_IG(end, :), t_final, constants.mu_EM);
rp_old = Xn_old_in(1:3) * constants.DU;
vp_old = Xn_old_in(4:6) * constants.VU;

% Build the active RTH target at the actual final epoch selected by the patch.
vInfTarget = local_build_seed_vinf_target_at_epoch(params, t_final);

% Sample the matching outgoing hyperbola and keep the periapsis closest to the current one.
[Rp_grid, Vp_grid] = rp_circle_from_vinfplus(vInfTarget, norm(rp_old), constants.muE, 720);
[~, best_idx] = min(vecnorm(Rp_grid - rp_old(:).', 2, 2));
rp_new = Rp_grid(best_idx, :);
vp_new = Vp_grid(best_idx, :);

% Preserve the incoming perigee velocity implied by the previous patch row.
vp_before = vp_old(:).' - FullSolution.deltavp(1, :);
FullSolution.deltavp(1, :) = vp_new - vp_before;
FullSolution.deltavpMag(1) = norm(FullSolution.deltavp(1, :));
if ismember('vpEarthFlyby', FullSolution.Properties.VariableNames)
    FullSolution.vpEarthFlyby(1, :) = vp_new;
end
if ismember('rpEarthFlyby', FullSolution.Properties.VariableNames)
    FullSolution.rpEarthFlyby(1, :) = rp_new;
end

% Store the corrected post-burn final state back in rotating coordinates.
nodes.X_nodes_IG(end, :) = in2rot([rp_new / constants.DU, vp_new / constants.VU], t_final, constants.mu_EM);
end

function needs_alignment = local_needs_free_final_rth_alignment(params)
%LOCAL_NEEDS_FREE_FINAL_RTH_ALIGNMENT Check whether the final-node seed must match an epoch-dependent RTH target.

needs_alignment = strcmpi(string(params.target_frame), "earth_sun_rth") && ...
    ~strcmpi(string(params.operational_mode), "fixed_arrival") && ...
    ~local_has_specified_final_time(params);
end

function t_final = local_free_final_rth_epoch_from_nodes(nodes, params)
%LOCAL_FREE_FINAL_RTH_EPOCH_FROM_NODES Read the final epoch as it will be represented by the packed NLP seed.

% Standard prescribed geometry pins t1 and keeps the patched segment durations.
if strcmpi(string(params.operational_mode), "standard") && strcmpi(string(params.geometry_mode), "prescribed_geometry")
    transfer_time = nodes.t_nodes_IG(end) - nodes.t_nodes_IG(1);
    t_final = params.initial_geometry_time + transfer_time;
    return
end

% Operational and free-geometry seeds already carry the absolute epoch.
t_final = nodes.t_nodes_IG(end);
end

function [n_nodes, info] = local_resolve_node_count(traj, params, constants)
%LOCAL_RESOLVE_NODE_COUNT Resolve the active multiple-shooting node count from the selected policy.

% Read the transfer duration in days from the assembled initial-guess timeline.
transfer_time_days = (traj.tTot(end) - traj.tTot(1)) * constants.TU / 86400;
node_count_mode = lower(string(params.node_count_mode));

% Use the adaptive node density when the requested policy is time-based.
if node_count_mode == "adaptive"
    n_nodes = ceil(params.nodes_per_day * transfer_time_days);
    n_nodes = max(params.n_nodes_min, n_nodes);
    n_nodes = min(params.n_nodes_max, n_nodes);
elseif node_count_mode == "adaptive_plus"
    [n_nodes, adaptive_plus_info] = local_resolve_adaptive_plus_count(traj, params, constants);
else
    n_nodes = params.n_nodes;
end

% Pack the node-count metadata used by logs, docs, and tests.
info = struct();
info.mode = char(node_count_mode);
info.transfer_time_days = transfer_time_days;
info.nodes_per_day = params.nodes_per_day;
info.n_nodes_base = n_nodes;
info.n_nodes_min = params.n_nodes_min;
info.n_nodes_max = params.n_nodes_max;
info.adaptive_plus_divisions = [];
info.adaptive_plus_upstream_nodes_requested = [];
info.adaptive_plus_upstream_nodes_used = [];
info.adaptive_plus_dt = [];
info.adaptive_plus_dt_days = [];
info.adaptive_plus_prefix_node_count = [];
info.adaptive_plus_even_start_idx_nodes = [];
info.adaptive_plus_between_internal_nodes = [];
if node_count_mode == "adaptive_plus"
    info.adaptive_plus_divisions = adaptive_plus_info.divisions;
    info.adaptive_plus_upstream_nodes_requested = adaptive_plus_info.upstream_nodes_requested;
    info.adaptive_plus_upstream_nodes_used = adaptive_plus_info.upstream_nodes_used;
    info.adaptive_plus_dt = adaptive_plus_info.dt;
    info.adaptive_plus_dt_days = adaptive_plus_info.dt_days;
    info.adaptive_plus_prefix_node_count = adaptive_plus_info.prefix_node_count;
    info.adaptive_plus_between_internal_nodes = adaptive_plus_info.between_internal_nodes;
end
end

function [n_prefix_nodes, info] = local_resolve_adaptive_plus_count(traj, params, constants)
%LOCAL_RESOLVE_ADAPTIVE_PLUS_COUNT Resolve the adaptive-plus spacing and the adaptive prefix node count.

% Read the two pericenter times used by the adaptive-plus construction.
t_start = traj.tTot(1);
t_moon_peri = traj.tTot(traj.moon_peri_idx);
t_earth_peri = traj.tTot(traj.earth_peri_idx);

% Normalize the user-selected adaptive-plus settings to integer counts.
divisions = max(1, round(params.adaptive_plus_divisions));
upstream_nodes_requested = max(0, round(params.adaptive_plus_upstream_nodes));

% Compute the uniform spacing imposed between Moon and Earth pericenters.
dt_adaptive_plus = (t_earth_peri - t_moon_peri) / divisions;
dt_adaptive_plus_days = abs(dt_adaptive_plus) * constants.TU / 86400;

% Clamp the upstream chain to the portion of the trajectory that exists before Moon peri.
upstream_nodes_used = upstream_nodes_requested;
if abs(dt_adaptive_plus) > eps
    max_upstream_nodes = floor(abs(t_moon_peri - t_start) / abs(dt_adaptive_plus));
    upstream_nodes_used = min(upstream_nodes_requested, max_upstream_nodes);
else
    upstream_nodes_used = 0;
end

% Identify the last evenly spaced node that remains in the adaptive prefix.
if upstream_nodes_used > 0
    t_prefix_end = t_moon_peri - upstream_nodes_used * dt_adaptive_plus;
else
    t_prefix_end = t_moon_peri;
end

% Compute the adaptive node count only on the remaining prefix duration.
prefix_duration_days = abs(t_prefix_end - t_start) * constants.TU / 86400;
if prefix_duration_days > 0
    n_prefix_nodes = ceil(params.nodes_per_day * prefix_duration_days);
    n_prefix_nodes = max(params.n_nodes_min, n_prefix_nodes);
    n_prefix_nodes = min(params.n_nodes_max, n_prefix_nodes);
else
    n_prefix_nodes = 1;
end

% Pack the adaptive-plus metadata reused by the node builder and tests.
info = struct();
info.divisions = divisions;
info.upstream_nodes_requested = upstream_nodes_requested;
info.upstream_nodes_used = upstream_nodes_used;
info.dt = dt_adaptive_plus;
info.dt_days = dt_adaptive_plus_days;
info.prefix_node_count = n_prefix_nodes;
info.between_internal_nodes = max(divisions - 1, 0);
end

function [t_nodes_IG, moon_peri_idx_nodes, earth_peri_idx_nodes, density_sep_idx] = local_build_standard_node_times(traj, params, n_nodes)
%LOCAL_BUILD_STANDARD_NODE_TIMES Build the legacy fixed/adaptive node times with the forced Moon-peri node.

% Reuse the legacy density split close to the end of the RPO leg.
RPO_idx = traj.RPO_idx;
density_sep_idx = round(RPO_idx * 0.99);

% Split the requested nodes into the nominal and dense portions.
dense_nodes = round(n_nodes * params.dense_nodes_multiplier);
nondense_nodes = n_nodes - dense_nodes;

% Keep at least two nodes in the nominal portion of the trajectory.
if nondense_nodes < 2
    nondense_nodes = max(nondense_nodes, 2);
    dense_nodes = max(n_nodes - nondense_nodes, 0);
end

% Build the initial node grid before the forced Moon-peri insertion.
t_nodes_IG_nondense = linspace(traj.tTot(1), traj.tTot(density_sep_idx), nondense_nodes).';
t_nodes_IG_dense = linspace(traj.tTot(density_sep_idx), traj.tTot(end), dense_nodes + 1).';
t_nodes_IG = [t_nodes_IG_nondense; t_nodes_IG_dense(2:end)];

% Insert the exact Moon-peri node into the standard node grid.
t_peri = traj.tTot(traj.moon_peri_idx);
idx_before = find(t_nodes_IG < t_peri, 1, 'last');
if isempty(idx_before)
    idx_before = 0;
end
moon_peri_idx_nodes = idx_before + 1;
earth_peri_idx_nodes = numel(t_nodes_IG) + 1;
t_nodes_IG = [t_nodes_IG(1:idx_before); t_peri; t_nodes_IG(idx_before+1:end)];
end

function [t_nodes_IG, moon_peri_idx_nodes, earth_peri_idx_nodes, density_sep_idx, info] = local_build_adaptive_plus_node_times(traj, params, constants, info)
%LOCAL_BUILD_ADAPTIVE_PLUS_NODE_TIMES Build the adaptive-plus node times with one evenly spaced Moon-to-Earth block.

% Read the pericenter times used to seed the evenly spaced late transfer nodes.
t_start = traj.tTot(1);
t_moon_peri = traj.tTot(traj.moon_peri_idx);
t_earth_peri = traj.tTot(traj.earth_peri_idx);
dt_adaptive_plus = info.adaptive_plus_dt;
upstream_nodes_used = info.adaptive_plus_upstream_nodes_used;

% Build the evenly spaced upstream chain that ends at Moon peri.
if upstream_nodes_used > 0
    t_upstream = (t_moon_peri - dt_adaptive_plus * (upstream_nodes_used:-1:1)).';
    t_prefix_end = t_upstream(1);
else
    t_upstream = [];
    t_prefix_end = t_moon_peri;
end

% Build the adaptive-plus prefix with the adaptive-derived count but without an extra dense tail at the handoff.
t_prefix_nodes = local_build_uniform_segment_times(t_start, t_prefix_end, info.adaptive_plus_prefix_node_count);
density_sep_idx = numel(t_prefix_nodes);

% Build the evenly spaced Moon-to-Earth internal nodes.
t_between_internal = (t_moon_peri + dt_adaptive_plus * (1:info.adaptive_plus_between_internal_nodes)).';

% Assemble the final adaptive-plus node timeline in chronological order.
if upstream_nodes_used > 0
    t_nodes_IG = [t_prefix_nodes; t_upstream(2:end); t_moon_peri; t_between_internal; t_earth_peri];
    info.adaptive_plus_even_start_idx_nodes = numel(t_prefix_nodes);
    moon_peri_idx_nodes = numel(t_prefix_nodes) + numel(t_upstream(2:end)) + 1;
else
    t_nodes_IG = [t_prefix_nodes; t_between_internal; t_earth_peri];
    info.adaptive_plus_even_start_idx_nodes = numel(t_prefix_nodes);
    moon_peri_idx_nodes = numel(t_prefix_nodes);
end

% Keep the Earth perigee at the last node of the adaptive-plus chain.
earth_peri_idx_nodes = numel(t_nodes_IG);
end

function t_nodes = local_build_dense_segment_times(t_start, t_end, n_nodes, dense_nodes_multiplier)
%LOCAL_BUILD_DENSE_SEGMENT_TIMES Build one adaptive-like dense segment between two selected times.

% Collapse degenerate segments to one endpoint node.
if abs(t_end - t_start) <= eps || n_nodes <= 1
    t_nodes = t_start;
    return
end

% Keep the two-point case on one simple uniform interval.
if n_nodes == 2
    t_nodes = linspace(t_start, t_end, 2).';
    return
end

% Split the segment into nominal and dense tails just like the adaptive mode.
dense_nodes = round(n_nodes * dense_nodes_multiplier);
nondense_nodes = n_nodes - dense_nodes;
if nondense_nodes < 2
    nondense_nodes = 2;
    dense_nodes = max(n_nodes - nondense_nodes, 0);
end

% Place the dense portion over the last one percent of the local segment.
t_dense_sep = t_start + 0.99 * (t_end - t_start);
t_nodes_nondense = linspace(t_start, t_dense_sep, nondense_nodes).';
t_nodes_dense = linspace(t_dense_sep, t_end, dense_nodes + 1).';
t_nodes = [t_nodes_nondense; t_nodes_dense(2:end)];
end

function t_nodes = local_build_uniform_segment_times(t_start, t_end, n_nodes)
%LOCAL_BUILD_UNIFORM_SEGMENT_TIMES Build one uniformly spaced segment between two selected times.

% Collapse degenerate segments to one endpoint node.
if abs(t_end - t_start) <= eps || n_nodes <= 1
    t_nodes = t_start;
    return
end

% Use one plain linear grid to avoid creating a dense cluster at the adaptive-plus handoff.
t_nodes = linspace(t_start, t_end, n_nodes).';
end

function impulses = init_impulses(FullSolution, nodes, constants, params)
%INIT_IMPULSES Construct the initial regularized impulse variables.
%
% PROTOTYPE:
%   impulses = init_impulses(FullSolution, nodes, constants, params)
%
% DESCRIPTION:
%   Map seed velocity jumps to U = (u, w, s) and record the alternate mapping
%   branches. Respect the active-node impulse policy and optional defect seed.
%
% INPUT:
%   FullSolution  Selected patch table with physical flyby states, phases,
%                 and times.
%   nodes  Shooting node states/times and their layout/encounter indices.
%   constants  Shared physical/scaling struct: DU [km], TU [s], VU [km/s],
%              mu_EM.
%   params  Transfer policy/settings struct; see the public RPO2DS_transfer
%           header.
%
% OUTPUT:
%   impulses  Regularized impulse seed and mapping-branch metadata.

U_nodes_IG = eps * ones(nodes.n_nodes, 3);
seed_impulse_mask = local_build_active_impulse_mask(nodes.n_nodes, params);
I_case2b = [];

% 1) Departure DV from RPO
DV0_CR3BP = FullSolution.dv0CR3BP(:);
[U0, is_case2b_0] = dv_to_regularized_U(DV0_CR3BP);
U_nodes_IG(1, :) = U0(:).';

if is_case2b_0
    I_case2b = [I_case2b; 1];
end

% Optionally initialize only interior impulses from CR3BP segment velocity defects.
if strcmpi(string(params.impulse_mode), "every_node") && local_use_defect_matching_impulse_seed(params)
    for j = 2:nodes.n_nodes - 1
        dt_seg = nodes.t_nodes_IG(j) - nodes.t_nodes_IG(j - 1);
        x_prev = nodes.X_nodes_IG(j - 1, :).';
        [x_prop, ~, ~, ~, ~] = propagate(x_prev, 0, dt_seg, constants.mu_EM, '3D', 'CR3BP', 'false');
        dv_node = nodes.X_nodes_IG(j, 4:6).' - x_prop(4:6);
        [U_node, is_case2b_node] = dv_to_regularized_U(dv_node);
        U_nodes_IG(j, :) = U_node(:).';
        if is_case2b_node
            I_case2b = [I_case2b; j]; %#ok<AGROW>
        end
    end
else
    seed_impulse_mask(2:end-1) = false;
end

% 2) Earth perigee DV from the patched-conic seed.
DVp_in = FullSolution.deltavp(:)';
xxp_in = rot2in(nodes.X_nodes_IG(end,:), nodes.t_nodes_IG(end), constants.mu_EM);
xxp_in = [xxp_in(1:3)*constants.DU, xxp_in(4:6)*constants.VU];
xxpbef_in = [xxp_in(1:3), xxp_in(4:6) - DVp_in];
xxpbef_rot = in2rot([xxpbef_in(1:3)/constants.DU, xxpbef_in(4:6)/constants.VU], nodes.t_nodes_IG(end), constants.mu_EM);
DVp_rot = nodes.X_nodes_IG(end,4:6) - xxpbef_rot(4:6);
[Up, is_case2b_p] = dv_to_regularized_U(DVp_rot);
U_nodes_IG(end, :) = Up(:).';
if is_case2b_p
    I_case2b = [I_case2b; nodes.n_nodes];
end

% Keep unseeded impulse nodes at the regularized zero-burn value.
U_nodes_IG(~seed_impulse_mask, :) = 0;
I_case2b = unique(I_case2b(:));

impulses = struct();
impulses.U_nodes_IG = U_nodes_IG;
impulses.I_case2b = I_case2b;
end

function use_seed = local_use_defect_matching_impulse_seed(params)
%LOCAL_USE_DEFECT_MATCHING_IMPULSE_SEED Read the optional high-feasibility impulse seeding policy.

use_seed = isfield(params, 'impulse_seed_mode') && ...
    strcmpi(string(params.impulse_seed_mode), "defect_matching");
end

function [y0, problem_data, meta, lb, ub] = pack_decision_vector(mode, vInfPlusEarth, nodes, impulses, constants, opts, params)
%PACK_DECISION_VECTOR Pack fixed-phase nodes, regularized impulses, and time variables.
%
% PROTOTYPE:
%   [y0, problem_data, meta, lb, ub] = pack_decision_vector(mode,
%   vInfPlusEarth, nodes, impulses, constants, opts, params)
%
% DESCRIPTION:
%   Build the shared node blocks and active time layout, then create variable
%   bounds and the model/target metadata consumed by the constraint routines.
%
% INPUT:
%   mode  opt minimizes regularized delta-v; reach maximizes v-infinity.
%   vInfPlusEarth  Earth-relative ECI target velocity/direction [km/s].
%   nodes  Shooting node states/times and their layout/encounter indices.
%   impulses  Initial regularized U variables and impulse-mapping branch
%             indices.
%   constants  Shared physical/scaling struct: DU [km], TU [s], VU [km/s],
%              mu_EM.
%   opts  Backend, iteration, tolerance, logging, plotting, and warm-start
%         options.
%   params  Transfer policy/settings struct; see the public RPO2DS_transfer
%           header.
%
% OUTPUT:
%   y0  Packed initial decision-vector column in the active node/time/phase
%       layout.
%   problem_data  Shared NLP constants, constraints policy, and
%                 packed-variable indices.
%   meta  Metadata required to decode the chosen packed layout after solving.
%   lb  Lower decision-variable bound column.
%   ub  Upper decision-variable bound column.

n_nodes = nodes.n_nodes;
time_layout = local_build_time_layout(nodes, constants.time_variable_mode, params);
active_impulse_mask = local_build_active_impulse_mask(n_nodes, params);

if mode == "opt"
    y0 = zeros(9 * n_nodes + time_layout.count, 1);
else
    y0 = zeros(9 * n_nodes + time_layout.count + 1, 1);
end

% Pack the state and regularized-impulse blocks on every node.
for j = 1:n_nodes
    offset = 9 * (j - 1);
    y0(offset + 1 : offset + 6) = nodes.X_nodes_IG(j, :).';
    y0(offset + 7 : offset + 9) = impulses.U_nodes_IG(j, :).';
end

% Pack either the absolute node times or the propagation durations after the node blocks.
y0(time_layout.time_idx) = time_layout.seed_values;

% Build the simple lower and upper bounds for any directly locked time variables.
[lb, ub] = local_build_time_bounds(numel(y0), time_layout);
[y0, lb, ub] = local_apply_impulse_policy_bounds(y0, lb, ub, active_impulse_mask);

% Build the shared problem payload used by objectives and constraints.
problem_data = build_base_problem_data(mode, vInfPlusEarth, nodes, impulses, constants, params, opts, time_layout, numel(y0));

% Store the setup metadata reused by post-processing and tests.
meta = struct();
meta = local_append_time_layout_meta(meta, problem_data);

if mode == "reach"
    vInfMag = norm(vInfPlusEarth);
    vinf_dir = vInfPlusEarth(:) / vInfMag;
    y0(end) = vInfMag;
    meta.vInfMag0 = vInfMag;
    meta.vinf_dir = vinf_dir;
end
end

function seed_time_info = build_seed_time_info(FullSolution, traj, params, guess_diagnostics, constants)
%BUILD_SEED_TIME_INFO Build one shared seed-time package for logs, tests, and result metadata.

seed_time_info = struct();
seed_time_info.mode = char(lower(string(params.time_mode)));
seed_time_info.seed_tof_lower_bound_days = NaN;
seed_time_info.seed_tof_exact_days = abs((traj.tTot(end) - traj.tTot(1)) * constants.TU / 86400);
seed_time_info.seed_tof_cap_days = params.maxTOF;
seed_time_info.seed_tof_keep = (seed_time_info.seed_tof_exact_days <= params.maxTOF);
seed_time_info.seed_flyby_inbound_days = NaN;
seed_time_info.max_tof_relax_factor = params.maxTOF_relax_factor;
seed_time_info.active_max_tof_days = params.maxTOF;
seed_time_info.max_tof_relaxed = false;
seed_time_info.summary_table = table();
seed_time_info.altitude_scan = table();

% Reuse the exact capped-seed metadata when the time-aware patch matcher produced it.
if ismember('seedToFLowerBoundDays', FullSolution.Properties.VariableNames)
    seed_time_info.seed_tof_lower_bound_days = FullSolution.seedToFLowerBoundDays(1);
end
if ismember('seedToFExactDays', FullSolution.Properties.VariableNames) && isfinite(FullSolution.seedToFExactDays(1))
    seed_time_info.seed_tof_exact_days = FullSolution.seedToFExactDays(1);
end
if ismember('seedToFCapDays', FullSolution.Properties.VariableNames) && isfinite(FullSolution.seedToFCapDays(1))
    seed_time_info.seed_tof_cap_days = FullSolution.seedToFCapDays(1);
end
if ismember('seedToFKeep', FullSolution.Properties.VariableNames)
    seed_time_info.seed_tof_keep = logical(FullSolution.seedToFKeep(1));
end
if ismember('seedFlybyInboundDays', FullSolution.Properties.VariableNames) && isfinite(FullSolution.seedFlybyInboundDays(1))
    seed_time_info.seed_flyby_inbound_days = FullSolution.seedFlybyInboundDays(1);
end

% Initial-state maxTOF is a total epoch-to-Earth window, not just transfer time.
if isfield(params, 'operational_mode') && strcmpi(string(params.operational_mode), "initial_state")
    seed_time_info.seed_tof_exact_days = abs((traj.tTot(end) - params.initial_epoch_time) * constants.TU / 86400);
    seed_time_info.seed_tof_cap_days = params.maxTOF;
    seed_time_info.seed_tof_keep = seed_time_info.seed_tof_exact_days <= params.maxTOF + 1e-12;
end

% Keep the altitude scan and refined screening summary for the driver and tests.
if isstruct(guess_diagnostics)
    if isfield(guess_diagnostics, 'best_patch_summary') && istable(guess_diagnostics.best_patch_summary)
        seed_time_info.summary_table = guess_diagnostics.best_patch_summary;
    end
    if isfield(guess_diagnostics, 'altitude_scan') && istable(guess_diagnostics.altitude_scan)
        seed_time_info.altitude_scan = guess_diagnostics.altitude_scan;
    end
end

% Relax the optimization TOF cap only when the selected max-time seed stays above the requested threshold.
if strcmpi(seed_time_info.mode, 'max_time') && ~seed_time_info.seed_tof_keep
    seed_time_info.active_max_tof_days = max(seed_time_info.seed_tof_cap_days, ...
        seed_time_info.seed_tof_exact_days * seed_time_info.max_tof_relax_factor);
    seed_time_info.max_tof_relaxed = seed_time_info.active_max_tof_days > seed_time_info.seed_tof_cap_days;
end
end

function phase_policy = local_update_operational_phase_seed(phase_policy, params, traj, constants)
%LOCAL_UPDATE_OPERATIONAL_PHASE_SEED Rebuild the fixed-arrival phase seed from the assembled transfer duration.

% Keep the other phase policies unchanged.
if ~isfield(params, 'operational_mode') || ~strcmpi(string(params.operational_mode), "fixed_arrival")
    return
end
if ~phase_policy.use_phase_variable
    return
end

% Use the assembled seed TOF so the packed equality matches the actual y0 timeline.
transfer_tof_nd = abs(traj.tTot(end) - traj.tTot(1));
phase_policy.phase_seed_nlp = params.initial_rpo_phase + ...
    (params.final_epoch_time - params.initial_epoch_time) - transfer_tof_nd;
end

function local_print_preoptimization_logs(FullSolution, guess_diagnostics, seed_time_info, traj, constants, opts)
%LOCAL_PRINT_PREOPTIMIZATION_LOGS Print the selected patch and initial-guess details before the optimization starts.

fprintf('\n--- Initial Guess Summary ---\n');
fprintf('%-32s %.6f days\n', 'Initial guess TOF:', abs((traj.tTot(end) - traj.tTot(1)) * constants.TU / 86400));
if ismember('deltavpMag', FullSolution.Properties.VariableNames)
    fprintf('%-32s %.6f km/s\n', 'Selected Earth-flyby dvp:', FullSolution.deltavpMag(1));
end
if ismember('T_lambert', FullSolution.Properties.VariableNames)
    fprintf('%-32s %.6f days\n', 'Selected Lambert TOF:', FullSolution.T_lambert(1) / 86400);
end
if ismember('T_RPOMoon', FullSolution.Properties.VariableNames)
    fprintf('%-32s %.6f days\n', 'Selected RPO-Moon TOF:', FullSolution.T_RPOMoon(1) * constants.TU / 86400);
end
fprintf('%-32s %.6f days\n', 'Requested max TOF:', seed_time_info.seed_tof_cap_days);
fprintf('%-32s %.6f days\n', 'Active optimization max TOF:', seed_time_info.active_max_tof_days);

if ~isstruct(guess_diagnostics) || ~isfield(guess_diagnostics, 'best_patch_diagnostics')
    return
end

patch_diag = guess_diagnostics.best_patch_diagnostics;
if isfield(patch_diag, 'filter_summary') && istable(patch_diag.filter_summary) && ~isempty(patch_diag.filter_summary)
    fprintf('\nPatch filter survival\n');
    disp(patch_diag.filter_summary);
end
if isfield(patch_diag, 'refinement_pool_summary') && istable(patch_diag.refinement_pool_summary) && ~isempty(patch_diag.refinement_pool_summary)
    fprintf('\nPatch refinement pool summary\n');
    disp(patch_diag.refinement_pool_summary);
end
if isfield(patch_diag, 'refinement_pool_table') && istable(patch_diag.refinement_pool_table) && ~isempty(patch_diag.refinement_pool_table)
    n_rows = min(local_patch_log_row_limit(opts), height(patch_diag.refinement_pool_table));
    fprintf('\nPatch refinement pool candidates (first %d rows)\n', n_rows);
    pool_columns = intersect({'candidateIndex','moonEarthIndex','rpoIndexLocal','rpoRowIndex','deltavpMag', ...
        'tofLowerBoundDays','tofExactDays','tofKeep','selected','poolRole'}, ...
        patch_diag.refinement_pool_table.Properties.VariableNames, 'stable');
    disp(patch_diag.refinement_pool_table(1:n_rows, pool_columns));
end
end

function n_rows = local_patch_log_row_limit(opts)
%LOCAL_PATCH_LOG_ROW_LIMIT Resolve the compact row cap for verbose patch candidate logs.

n_rows = opts.patch_log_candidate_rows;
if ~isfinite(n_rows)
    n_rows = inf;
else
    n_rows = max(1, round(n_rows));
end
end

function phase_policy = finalize_phase_policy(phase_policy, FullSolution)
%FINALIZE_PHASE_POLICY Finalize the phase metadata after the patched seed row has been selected.

% Read the selected seed phase from the patched solution table.
phase_seed = 0;
    if RPO2DS_has_patch_solution(FullSolution) && ismember('phase', FullSolution.Properties.VariableNames)
        phase_seed = FullSolution.phase(1);
    end

% Keep the requested and seeded phase values together in the shared metadata.
phase_policy.phase_seed = phase_seed;
phase_policy.phase_seed_wrapped = phase_policy.phase_origin + mod(phase_seed - phase_policy.phase_origin, phase_policy.phase_period);
end

function [y0, problem_data, meta, lb, ub] = build_setup_data(mode, vInfPlusEarth, FullSolution, nodes, impulses, constants, params, opts, phase_policy, seed_time_info)
%BUILD_SETUP_DATA Build the active phase-policy NLP setup and reusable metadata.
%
% PROTOTYPE:
%   [y0, problem_data, meta, lb, ub] = build_setup_data(mode, vInfPlusEarth,
%   FullSolution, nodes, impulses, constants, params, opts, phase_policy,
%   seed_time_info)
%
% DESCRIPTION:
%   Select fixed/phase-aware packing, apply the Sun-angle and time policies,
%   and attach the seed/operational metadata used by solving and decoding.
%
% INPUT:
%   mode  opt minimizes regularized delta-v; reach maximizes v-infinity.
%   vInfPlusEarth  Earth-relative ECI target velocity/direction [km/s].
%   FullSolution  Selected patch table with physical flyby states, phases,
%                 and times.
%   nodes  Shooting node states/times and their layout/encounter indices.
%   impulses  Initial regularized U variables and impulse-mapping branch
%             indices.
%   constants  Shared physical/scaling struct: DU [km], TU [s], VU [km/s],
%              mu_EM.
%   params  Transfer policy/settings struct; see the public RPO2DS_transfer
%           header.
%   opts  Backend, iteration, tolerance, logging, plotting, and warm-start
%         options.
%   phase_policy  Database phase grid/origin/period [TU] and active
%                 phase-source policy.
%   seed_time_info  Seed flight-time bounds, selected cap, and exact timing
%                   metadata.
%
% OUTPUT:
%   y0  Packed initial decision-vector column in the active node/time/phase
%       layout.
%   problem_data  Shared NLP constants, constraints policy, and
%                 packed-variable indices.
%   meta  Metadata required to decode the chosen packed layout after solving.
%   lb  Lower decision-variable bound column.
%   ub  Upper decision-variable bound column.

if ~phase_policy.use_phase_variable
    [y0, problem_data, meta, lb, ub] = pack_decision_vector(mode, vInfPlusEarth, nodes, impulses, constants, opts, params);
    [y0, problem_data, meta, lb, ub] = local_apply_theta_sun0_variable(mode, y0, problem_data, meta, lb, ub, params);
    meta.phase_mode = char(phase_policy.mode);
    meta.use_phase_variable = false;
    meta.phase_seed = phase_policy.phase_seed;
    meta.phase_request = phase_policy.phase_request_wrapped;
    meta.phase_shortlist_rows = phase_policy.shortlist_rows;
    meta.phase_shortlist_phases = phase_policy.shortlist_phases;
    meta.phase_shortlist_distances = phase_policy.shortlist_distances;
    meta = local_append_time_layout_meta(meta, problem_data);
    meta.time_mode = seed_time_info.mode;
    meta.seed_tof_lower_bound_days = seed_time_info.seed_tof_lower_bound_days;
    meta.seed_tof_exact_days = seed_time_info.seed_tof_exact_days;
    meta.seed_tof_cap_days = seed_time_info.seed_tof_cap_days;
    meta.seed_tof_keep = seed_time_info.seed_tof_keep;
    meta.seed_flyby_inbound_days = seed_time_info.seed_flyby_inbound_days;
    meta.max_tof_relax_factor = seed_time_info.max_tof_relax_factor;
    meta.active_max_tof_days = seed_time_info.active_max_tof_days;
    meta.max_tof_relaxed = seed_time_info.max_tof_relaxed;
    meta.decision_vector_length = numel(y0);
    problem_data.maxTOF = seed_time_info.active_max_tof_days * 86400 / constants.TU;
    problem_data.maxTransferTOF = problem_data.maxTOF;
    return
end

% Build the new phase-aware decision vector and its NLP/problem structure.
[y0, base_problem_data, lb, ub, phase_theta_idx] = pack_decision_vector_phase(mode, vInfPlusEarth, FullSolution, nodes, impulses, constants, params, opts, phase_policy);
problem_data = build_phase_problem_data(base_problem_data, phase_policy, phase_theta_idx, numel(y0));
meta = struct();
[y0, problem_data, meta, lb, ub] = local_apply_theta_sun0_variable(mode, y0, problem_data, meta, lb, ub, params);

% Pack the phase-aware metadata reused by post-processing and tests.
meta.phase_mode = char(phase_policy.mode);
meta.use_phase_variable = true;
meta.phase_source = problem_data.phase_source;
meta.has_phase_variable = problem_data.has_phase_variable;
meta.phase_seed = phase_policy.phase_seed;
meta.phase_request = phase_policy.phase_request_wrapped;
meta.phase_origin = phase_policy.phase_origin;
meta.phase_period = phase_policy.phase_period;
meta.phase_shortlist_rows = phase_policy.shortlist_rows;
meta.phase_shortlist_phases = phase_policy.shortlist_phases;
meta.phase_shortlist_distances = phase_policy.shortlist_distances;
meta.phase_theta_idx = phase_theta_idx;
if ~isempty(phase_theta_idx)
    meta.phase_bounds = [lb(phase_theta_idx), ub(phase_theta_idx)];
else
    meta.phase_bounds = [];
end
meta = local_append_time_layout_meta(meta, problem_data);
meta.time_mode = seed_time_info.mode;
meta.seed_tof_lower_bound_days = seed_time_info.seed_tof_lower_bound_days;
meta.seed_tof_exact_days = seed_time_info.seed_tof_exact_days;
meta.seed_tof_cap_days = seed_time_info.seed_tof_cap_days;
meta.seed_tof_keep = seed_time_info.seed_tof_keep;
meta.seed_flyby_inbound_days = seed_time_info.seed_flyby_inbound_days;
meta.max_tof_relax_factor = seed_time_info.max_tof_relax_factor;
meta.active_max_tof_days = seed_time_info.active_max_tof_days;
meta.max_tof_relaxed = seed_time_info.max_tof_relaxed;
meta.decision_vector_length = numel(y0);
if strcmp(params.operational_mode, 'standard')
    problem_data.maxTOF = seed_time_info.active_max_tof_days * 86400 / constants.TU;
    problem_data.maxTransferTOF = problem_data.maxTOF;
elseif strcmp(params.operational_mode, 'initial_state')
    problem_data.maxTotalTime = seed_time_info.active_max_tof_days * 86400 / constants.TU;
end
if mode == "reach"
    meta.vinf_mag_idx = numel(y0);
end
end

function [y0, problem_data, meta, lb, ub] = local_apply_theta_sun0_variable(mode, y0, problem_data, meta, lb, ub, params)
%LOCAL_APPLY_THETA_SUN0_VARIABLE Add thetaSun0 to the NLP only when the target frame requests it.

% Store the target-frame controls in every problem payload.
problem_data.target_frame = char(lower(string(params.target_frame)));
problem_data.theta_sun0_mode = char(lower(string(params.theta_sun0_mode)));
problem_data.theta_sun0_value = params.theta_sun0_value;
problem_data.theta_sun0_seed = params.theta_sun0_seed;
problem_data.theta_sun0_idx = [];
problem_data.ws_sun_em_nd = params.ws_sun_em_nd;
problem_data.vinf_target_rth = params.vinf_target_rth_km_s(:);
problem_data.T_RPO = -pi / params.ws_sun_em_nd;
problem_data.initial_wait_cap = 2 * problem_data.T_RPO;

% Store matching metadata for result diagnostics and tests.
meta.target_frame = problem_data.target_frame;
meta.theta_sun0_mode = problem_data.theta_sun0_mode;
meta.theta_sun0_value = problem_data.theta_sun0_value;
meta.theta_sun0_seed = problem_data.theta_sun0_seed;
meta.theta_sun0_idx = [];
meta.ws_sun_em_nd = problem_data.ws_sun_em_nd;
meta.vinf_target_rth = problem_data.vinf_target_rth;
meta.T_RPO = problem_data.T_RPO;
meta.initial_wait_cap = problem_data.initial_wait_cap;

% Keep thetaSun0 as a parameter unless the RTH target needs a free Sun phase.
if lower(string(params.theta_sun0_mode)) ~= "free" || lower(string(params.target_frame)) ~= "earth_sun_rth"
    problem_data.nvars = numel(y0);
    return
end

% Insert thetaSun0 before vInfMag in reach mode and at the end in opt mode.
if mode == "reach"
    theta_idx = numel(y0);
    y0 = [y0(1:end-1); params.theta_sun0_seed; y0(end)];
    lb = [lb(1:end-1); 0; lb(end)];
    ub = [ub(1:end-1); 2 * pi; ub(end)];
    problem_data.vinf_mag_idx = numel(y0);
else
    y0 = [y0; params.theta_sun0_seed];
    lb = [lb; 0];
    ub = [ub; 2 * pi];
    theta_idx = numel(y0);
end

% Expose the packed thetaSun0 index to dense, CasADi, and reports.
problem_data.theta_sun0_idx = theta_idx;
problem_data.nvars = numel(y0);
meta.theta_sun0_idx = theta_idx;
meta.theta_sun0_bounds = [lb(theta_idx), ub(theta_idx)];
end

function [y0, base_problem_data, lb, ub, phase_theta_idx] = pack_decision_vector_phase(mode, vInfPlusEarth, ~, nodes, impulses, constants, params, opts, phase_policy)
%PACK_DECISION_VECTOR_PHASE Pack the phase-aware node/time/phase decision-vector layout.
%
% PROTOTYPE:
%   [y0, base_problem_data, lb, ub, phase_theta_idx] = ...
%       pack_decision_vector_phase(mode, vInfPlusEarth, ~, nodes, impulses, ...
%       constants, params, opts, phase_policy)
%
% DESCRIPTION:
%   Keep fixed-phase node blocks [X(6); U(3)]. Add an explicit phase variable
%   only when its source requires one; time-derived phase uses epoch fields.
%
% INPUT:
%   mode  opt minimizes regularized delta-v; reach maximizes v-infinity.
%   vInfPlusEarth  Earth-relative ECI target velocity/direction [km/s].
%   ~  Unused argument retained by the shared calling interface.
%   nodes  Shooting state/time arrays, layout, and encounter-node indices.
%   impulses  Regularized U seed and impulse-mapping branch indices.
%   constants  Physical/scaling constants: DU [km], TU [s], VU [km/s], mu_EM.
%   params  Normalized transfer policies; see the public RPO2DS_transfer
%           header.
%   opts  Backend, iteration, tolerance, display, and warm-start options.
%   phase_policy  RPO database phase grid/period [TU] and the active phase
%                 source.
%
% OUTPUT:
%   y0  Packed initial decision-vector column.
%   base_problem_data  Base model, target, and time-layout struct.
%   lb  Lower decision-variable bound column.
%   ub  Upper decision-variable bound column.
%   phase_theta_idx  Explicit phase-variable index, when that layout uses
%                    one.

[y_legacy, base_problem_data, meta_legacy, lb_legacy, ub_legacy] = pack_decision_vector(mode, vInfPlusEarth, nodes, impulses, constants, opts, params);
phase_seed = phase_policy.phase_seed_wrapped;
if isfield(phase_policy, 'phase_seed_nlp') && ~isempty(phase_policy.phase_seed_nlp)
    phase_seed = phase_policy.phase_seed_nlp;
end
phase_variable_mode = "explicit";
if isfield(phase_policy, 'phase_variable_mode') && ~isempty(phase_policy.phase_variable_mode)
    phase_variable_mode = string(phase_policy.phase_variable_mode);
end

% Initial-state mode derives the RPO phase from t(1), so no theta variable is packed.
if phase_variable_mode == "from_initial_time"
    y0 = y_legacy;
    lb = lb_legacy;
    ub = ub_legacy;
    phase_theta_idx = [];
    phase_seed_geometry = local_wrap_phase_value(phase_seed, phase_policy.phase_origin, phase_policy.phase_period);
    y0(1:6) = local_build_phase_consistent_first_node_from_theta(phase_seed_geometry, y0, constants.xx0RPOBase(:), impulses.I_case2b(:), constants.mu_EM);
    return
end

% Append theta_rpo before vInfMag in reachability mode and at the end in opt mode.
if mode == "opt"
    y0 = [y_legacy; phase_seed];
    phase_theta_idx = numel(y0);
else
    y0 = [y_legacy(1:end-1); phase_seed; y_legacy(end)];
    phase_theta_idx = numel(y0) - 1;
end

% Build default bounds and clamp the phase when the user requested other_phase.
lb = local_insert_phase_variable(lb_legacy, phase_theta_idx, mode, -inf);
ub = local_insert_phase_variable(ub_legacy, phase_theta_idx, mode, inf);
if isfield(phase_policy, 'phase_bounds_nlp') && ~isempty(phase_policy.phase_bounds_nlp)
    lb(phase_theta_idx) = phase_policy.phase_bounds_nlp(1);
    ub(phase_theta_idx) = phase_policy.phase_bounds_nlp(2);
else
    lb(phase_theta_idx) = phase_policy.phase_origin;
    ub(phase_theta_idx) = phase_policy.phase_origin + phase_policy.phase_period;
end
if phase_policy.mode == "other_phase"
    y0(phase_theta_idx) = phase_policy.phase_request_wrapped;
    lb(phase_theta_idx) = phase_policy.phase_request_wrapped;
    ub(phase_theta_idx) = phase_policy.phase_request_wrapped;
end

% Keep the reachability vInf seed unchanged after inserting theta_rpo.
if mode == "reach"
    y0(end) = meta_legacy.vInfMag0;
end

% Rebuild the first node so it is consistent with the active departure phase.
y0(1:6) = local_build_phase_consistent_first_node(y0, phase_theta_idx, constants.xx0RPOBase(:), impulses.I_case2b(:), constants.mu_EM);
end

function problem_data = build_base_problem_data(mode, vInfPlusEarth, nodes, impulses, constants, params, opts, time_layout, nvars)
%BUILD_BASE_PROBLEM_DATA Assemble shared model, target, and time-layout constraint data.
%
% PROTOTYPE:
%   problem_data = build_base_problem_data(mode, vInfPlusEarth, nodes, ...
%       impulses, constants, params, opts, time_layout, nvars)
%
% DESCRIPTION:
%   Record physical scales, normalized primary/time limits, periapsis node,
%   impulse branches, target frame, and packed-variable indices.
%
% INPUT:
%   mode / vInfPlusEarth       opt/reach and Earth-relative ECI target [km/s].
%   nodes / impulses           Shooting states/times and regularized U seed.
%   constants / params / opts  Physical scales and normalized solve policies.
%   time_layout                Packed time indices, anchors, and locks [TU].
%   nvars                      Total decision-vector length.
%
% OUTPUT:
%   problem_data  Shared NLP model/target/constants and variable indices.

problem_data = struct();
problem_data.mode = char(mode);
problem_data.mu = constants.mu_EM;
problem_data.rminEarth = constants.rminEarth_DU;
problem_data.rminMoon = constants.rminMoon_DU;
problem_data.maxFinalEarthDistance = [];
if ~isempty(params.max_final_earth_distance_km)
    problem_data.maxFinalEarthDistance = params.max_final_earth_distance_km / constants.DU;
end
problem_data.muE = constants.muE;
problem_data.DU = constants.DU;
problem_data.VU = constants.VU;
problem_data.TU = constants.TU;
problem_data.n = nodes.n_nodes;
problem_data.maxTOF = constants.maxTOF;
problem_data.idx_peri_moon = nodes.moon_peri_idx_nodes;
problem_data.I_2b = impulses.I_case2b(:);
problem_data.impulse_mode = char(lower(string(params.impulse_mode)));
problem_data.active_impulse_mask = local_build_active_impulse_mask(nodes.n_nodes, params);
problem_data.x0_rpo_base = constants.xx0RPOBase(:);
problem_data.x_ini_fixed = constants.xx0RPO(:);
    problem_data.use_phase_variable = false;
    problem_data.time_variable_mode = time_layout.mode; % Store the active time layout used by the NLP
    problem_data.time_anchor_mode = time_layout.time_anchor_mode; % Store whether propagation times are anchored at t1 or tN
    problem_data.geometry_mode = time_layout.geometry_mode;
    problem_data.t0_anchor = time_layout.t0_anchor; % Keep the fixed anchor time used by propagation_times
    problem_data.time_idx = time_layout.time_idx(:); % Store the packed time-variable indices inside y
    problem_data.dt_idx = time_layout.dt_idx(:);
    problem_data.t1_idx = time_layout.t1_idx;
problem_data.time_count = time_layout.count; % Keep the number of packed time variables for this layout
problem_data.final_time_idx = time_layout.final_time_idx; % Track the final packed time-variable index
problem_data.initial_geometry_time = time_layout.initial_geometry_time;
problem_data.lock_t_node1 = time_layout.lock_t_node1;
problem_data.lock_t_nodeN = time_layout.lock_t_nodeN;
problem_data.t_node1_fixed = time_layout.t_node1_fixed;
problem_data.t_nodeN_fixed = time_layout.t_nodeN_fixed;
problem_data.has_t1_variable = time_layout.has_t1_variable;
problem_data.operational_mode = char(params.operational_mode);
problem_data.initial_epoch_time = params.initial_epoch_time;
problem_data.initial_rpo_phase = params.initial_rpo_phase;
problem_data.final_epoch_time = params.final_epoch_time;
problem_data.maxTransferTOF = constants.maxTOF;
problem_data.maxTotalTime = [];
problem_data.available_time_nd = [];
if ~isempty(params.final_epoch_time)
    problem_data.available_time_nd = params.final_epoch_time - params.initial_epoch_time;
end
if strcmp(params.operational_mode, 'initial_state')
    problem_data.maxTransferTOF = params.maxTransferTOF * 86400 / constants.TU;
    problem_data.maxTOF = problem_data.maxTransferTOF;
    problem_data.maxTotalTime = params.maxTOF * 86400 / constants.TU;
elseif strcmp(params.operational_mode, 'fixed_arrival')
    problem_data.maxTransferTOF = params.maxTransferTOF * 86400 / constants.TU;
    problem_data.maxTOF = problem_data.maxTransferTOF;
end
if strcmp(params.operational_mode, 'fixed_arrival') && ~isempty(problem_data.available_time_nd)
    problem_data.maxTransferTOF = min(problem_data.maxTransferTOF, problem_data.available_time_nd);
    problem_data.maxTOF = problem_data.maxTransferTOF;
end
    problem_data.nvars = nvars;
if mode == "opt"
    problem_data.vinf_target = vInfPlusEarth(:);
else
    problem_data.vinf_dir = vInfPlusEarth(:) / norm(vInfPlusEarth);
    problem_data.dv_budget = opts.dv_budget;
    problem_data.vinf_mag_idx = nvars;
end
end

function problem_data = build_phase_problem_data(base_problem_data, phase_policy, phase_theta_idx, nvars)
%BUILD_PHASE_PROBLEM_DATA Attach departure-phase constraints to the shared problem data.
%
% PROTOTYPE:
%   problem_data = build_phase_problem_data(base_problem_data, ...
%       phase_policy, phase_theta_idx, nvars)
%
% DESCRIPTION:
%   Retain the base constraint schema and add the RPO departure-state/phase
%   origin, period, source, and bounds used by the phase-aware boundary.
%
% INPUT:
%   base_problem_data  Fixed-phase model/target/time schema.
%   phase_policy       RPO grid/period and active phase-source policy [TU].
%   phase_theta_idx    Index of the explicit phase variable, when used.
%   nvars              Total packed decision-vector length.
%
% OUTPUT:
%   problem_data       Phase-aware NLP model and decoding fields.

problem_data = base_problem_data;
problem_data.use_phase_variable = true;
problem_data.has_phase_variable = ~isempty(phase_theta_idx);
problem_data.phase_source = 'variable';
if isfield(phase_policy, 'phase_variable_mode') && strcmp(string(phase_policy.phase_variable_mode), "from_initial_time")
    problem_data.phase_source = 'initial_time';
end
problem_data.phase_theta_idx = phase_theta_idx;
problem_data.phase_origin = phase_policy.phase_origin;
problem_data.phase_period = phase_policy.phase_period;
problem_data.phase_mode = char(phase_policy.mode);
problem_data.phase_seed = phase_policy.phase_seed_wrapped;
if isfield(phase_policy, 'phase_seed_nlp') && ~isempty(phase_policy.phase_seed_nlp)
    problem_data.phase_seed = phase_policy.phase_seed_nlp;
end
problem_data.nvars = nvars;
if strcmp(problem_data.mode, 'reach')
    problem_data.vinf_mag_idx = nvars;
end
end

function time_layout = local_build_time_layout(nodes, time_variable_mode, params)
%LOCAL_BUILD_TIME_LAYOUT Define absolute-node or propagation-duration time variables.
%
% PROTOTYPE:
%   time_layout = local_build_time_layout(nodes, time_variable_mode, params)
%
% DESCRIPTION:
%   Choose the initial/final time anchor, active epoch variable, and fixed-time
%   locks. Store indices needed to reconstruct absolute epochs in TU.
%
% INPUT:
%   nodes               Shooting node count and seed epochs [TU].
%   time_variable_mode  node_times or propagation_times.
%   params              Operational/geometry policies and prescribed epochs.
%
% OUTPUT:
%   time_layout         Struct with time/duration/epoch indices, anchors,
%                        locks, fixed values, and number of time variables.

time_layout = struct();
time_layout.mode = char(lower(string(time_variable_mode)));
time_layout.time_anchor_mode = 'initial';
time_layout.geometry_mode = 'free_geometry';
time_layout.initial_geometry_time = [];
time_layout.t_node1_fixed = [];
time_layout.t_nodeN_fixed = [];
time_layout.lock_t_node1 = false;
time_layout.lock_t_nodeN = false;
time_layout.has_t1_variable = false;
time_layout.t1_idx_local = [];
time_layout.dt_idx_local = [];
time_seed = nodes.t_nodes_IG(:);
dt_seed = diff(time_seed);
time_layout.t0_anchor = time_seed(1);

% Keep the prescribed-geometry lock only in the standard departure-driven mode.
if strcmpi(string(params.operational_mode), "standard")
    time_layout.geometry_mode = char(lower(string(params.geometry_mode)));
    time_layout.initial_geometry_time = params.initial_geometry_time;
end

% Lock the departure epoch when the geometry must stay aligned to one selected phase.
if strcmp(time_layout.geometry_mode, 'prescribed_geometry')
    time_layout.lock_t_node1 = true;
    time_layout.t_node1_fixed = params.initial_geometry_time;
end

% Lock the Earth-perigee epoch when the target frame uses a specified final time.
if local_has_specified_final_time(params)
    time_layout.lock_t_nodeN = true;
    time_layout.t_nodeN_fixed = params.final_time_seed;
end

if strcmp(time_layout.mode, 'propagation_times')
    if strcmpi(string(params.operational_mode), "fixed_arrival")
        % Anchor fixed-arrival propagation times at the prescribed final epoch.
        time_layout.time_anchor_mode = 'final';
        time_layout.lock_t_nodeN = true;
        time_layout.t_nodeN_fixed = params.final_epoch_time;
        time_layout.t0_anchor = params.final_epoch_time;
        time_layout.seed_values = dt_seed;
    else
        % Keep the propagation-times layout and seed the enforced departure epoch directly.
        time_layout.has_t1_variable = true;
        t1_seed = time_seed(1);
        if time_layout.lock_t_node1
            t1_seed = time_layout.t_node1_fixed;
        end
        time_layout.t0_anchor = t1_seed;
        time_layout.seed_values = [t1_seed; dt_seed];
    end
else
    time_layout.mode = 'node_times';
    % Shift the whole seed timeline so the first node already matches the prescribed epoch.
    if time_layout.lock_t_node1
        time_seed = time_seed - time_seed(1) + time_layout.t_node1_fixed;
    end
    time_layout.seed_values = time_seed; % Keep the absolute node-time layout with the enforced first epoch
    time_layout.has_t1_variable = true;
    time_layout.t0_anchor = time_seed(1);
end
time_layout.count = numel(time_layout.seed_values);
time_layout.time_idx = (9 * nodes.n_nodes + 1 : 9 * nodes.n_nodes + time_layout.count).';
if strcmp(time_layout.mode, 'propagation_times')
    if time_layout.has_t1_variable
        time_layout.t1_idx_local = 1;
        if time_layout.count > 1
            time_layout.dt_idx_local = 2:time_layout.count;
        end
    else
        time_layout.dt_idx_local = 1:time_layout.count;
    end
end
if isempty(time_layout.dt_idx_local)
    time_layout.dt_idx = [];
else
    time_layout.dt_idx = time_layout.time_idx(time_layout.dt_idx_local);
end
if isempty(time_layout.t1_idx_local)
    time_layout.t1_idx = [];
else
    time_layout.t1_idx = time_layout.time_idx(time_layout.t1_idx_local);
end
if strcmp(time_layout.mode, 'node_times')
    time_layout.final_time_idx = time_layout.time_idx(end);
else
    time_layout.final_time_idx = [];
end
end

function [lb, ub] = local_build_time_bounds(nvars, time_layout)
%LOCAL_BUILD_TIME_BOUNDS Build the time bounds, including the enforced departure epoch when requested.

lb = -inf(nvars, 1);
ub = inf(nvars, 1);

% Pin the first node time when the prescribed geometry requires one fixed initial epoch.
if time_layout.lock_t_node1
    if strcmp(time_layout.mode, 'propagation_times')
        lb(time_layout.t1_idx) = time_layout.t_node1_fixed;
        ub(time_layout.t1_idx) = time_layout.t_node1_fixed;
    else
        lb(time_layout.time_idx(1)) = time_layout.t_node1_fixed;
        ub(time_layout.time_idx(1)) = time_layout.t_node1_fixed;
    end
end

% Pin the final absolute node time when it is stored directly in node_times mode.
if time_layout.lock_t_nodeN && strcmp(time_layout.mode, 'node_times')
    lb(time_layout.time_idx(end)) = time_layout.t_nodeN_fixed;
    ub(time_layout.time_idx(end)) = time_layout.t_nodeN_fixed;
end

end

function active_impulse_mask = local_build_active_impulse_mask(n_nodes, params)
%LOCAL_BUILD_ACTIVE_IMPULSE_MASK Build the node-wise impulse activation mask from the selected policy.

active_impulse_mask = true(n_nodes, 1);
if strcmpi(string(params.impulse_mode), "endpoints_only")
    active_impulse_mask(:) = false;
    active_impulse_mask([1, n_nodes]) = true;
end
end

function [y0, lb, ub] = local_apply_impulse_policy_bounds(y0, lb, ub, active_impulse_mask)
%LOCAL_APPLY_IMPULSE_POLICY_BOUNDS Pin inactive regularized impulse variables to zero.

for j = 1:numel(active_impulse_mask)
    if active_impulse_mask(j)
        continue
    end
    impulse_idx = 9 * (j - 1) + (7:9);
    y0(impulse_idx) = 0;
    lb(impulse_idx) = 0;
    ub(impulse_idx) = 0;
end
end

function bounds_out = local_insert_phase_variable(bounds_in, phase_theta_idx, mode, inserted_value)
%LOCAL_INSERT_PHASE_VARIABLE Insert the phase variable bound in the same position used by theta_rpo.

if mode == "opt"
    bounds_out = [bounds_in; inserted_value];
else
    bounds_out = [bounds_in(1:end-1); inserted_value; bounds_in(end)];
end

% Keep the inserted scalar at the expected phase index for clarity.
bounds_out(phase_theta_idx) = inserted_value;
end

function meta = local_append_time_layout_meta(meta, problem_data)
%LOCAL_APPEND_TIME_LAYOUT_META Copy the shared time-layout metadata into the setup/result meta struct.

meta.time_variable_mode = problem_data.time_variable_mode;
meta.time_anchor_mode = problem_data.time_anchor_mode;
meta.geometry_mode = problem_data.geometry_mode;
meta.t0_anchor = problem_data.t0_anchor;
meta.time_idx = problem_data.time_idx(:);
meta.dt_idx = problem_data.dt_idx(:);
meta.t1_idx = problem_data.t1_idx;
meta.time_block_count = problem_data.time_count;
meta.final_time_idx = problem_data.final_time_idx;
meta.initial_geometry_time = problem_data.initial_geometry_time;
meta.lock_t_node1 = problem_data.lock_t_node1;
meta.lock_t_nodeN = problem_data.lock_t_nodeN;
meta.t_node1_fixed = problem_data.t_node1_fixed;
meta.t_nodeN_fixed = problem_data.t_nodeN_fixed;
meta.impulse_mode = problem_data.impulse_mode;
meta.active_impulse_mask = problem_data.active_impulse_mask;
end

function meta = local_append_operational_meta(meta, problem_data, operational_guess)
%LOCAL_APPEND_OPERATIONAL_META Copy the operational-mode metadata into the setup/result meta struct.

meta.operational_mode = problem_data.operational_mode;
meta.initial_epoch_time = problem_data.initial_epoch_time;
meta.initial_rpo_phase = problem_data.initial_rpo_phase;
meta.final_epoch_time = problem_data.final_epoch_time;
meta.available_time_nd = problem_data.available_time_nd;
meta.seed_candidate_index = [];
meta.operational_screening_summary = table();
if isstruct(operational_guess) && isfield(operational_guess, 'selected_row')
    meta.seed_candidate_index = operational_guess.selected_row;
    meta.operational_screening_summary = operational_guess.screening;
end
end

function [y0, problem_data, meta, lb, ub, setup_data] = local_build_continuation_setup(base_result, vInfPlusEarth, params, opts)
%LOCAL_BUILD_CONTINUATION_SETUP Rebuild a reusable NLP setup from a previous transfer result.
%
% PROTOTYPE:
%   [y0, problem_data, meta, lb, ub, setup_data] = ...
%       local_build_continuation_setup(base_result, vInfPlusEarth, params, opts)
%
% DESCRIPTION:
%   Retain node/impulse structure and the previous solution seed while updating
%   the requested target and initial geometry-time lock.
%
% INPUT:
%   base_result    Previous RPO2DS_transfer result with reusable setup data.
%   vInfPlusEarth  Updated Earth-relative ECI target [km/s].
%   params / opts Normalized policies and continuation/warm-start settings.
%
% OUTPUT:
%   y0, lb, ub    Packed seed and variable bounds in their original units.
%   problem_data  Updated NLP model/target/time/phase fields.
%   meta          Updated decoding metadata.
%   setup_data    Reusable seed, bound, problem, and metadata package.

if ~isfield(base_result, 'y_opt') || isempty(base_result.y_opt)
    error('RPO2DS_transfer:MissingContinuationYOpt', ...
        'opts.continuation_base_result must contain one previous optimized decision vector.');
end

% Start from the previous optimized decision vector unless the caller overrides it explicitly.
y0 = base_result.y_opt(:);
if isfield(opts, 'y0_override') && ~isempty(opts.y0_override)
    y0 = opts.y0_override(:);
end

% Reuse the previous problem layout, metadata, and bounds directly.
problem_data = base_result.setup.problem_data;
meta = base_result.meta;
lb = base_result.setup.lb(:);
ub = base_result.setup.ub(:);

% Update the enforced initial-time lock only for prescribed-geometry continuation.
if strcmpi(string(params.geometry_mode), "prescribed_geometry")
    [y0, problem_data, meta, lb, ub] = local_update_continuation_time_lock(y0, problem_data, meta, lb, ub, params.initial_geometry_time);
end

% Update the target vector, reachability direction, and budget for direction sweeps.
[y0, problem_data, meta] = local_update_continuation_target(y0, problem_data, meta, vInfPlusEarth, params, opts);

% Pack the continuation setup used by plots, tests, and diagnostics.
setup_data = struct();
setup_data.y0_seed = y0;
setup_data.y0 = y0;
setup_data.problem_data = problem_data;
setup_data.lb = lb;
setup_data.ub = ub;
end

function [y0, problem_data, meta] = local_update_continuation_target(y0, problem_data, meta, vInfPlusEarth, params, opts)
%LOCAL_UPDATE_CONTINUATION_TARGET Update the fixed target or reachability direction while preserving the NLP layout.

% Keep the target-frame metadata aligned with the current call settings.
problem_data.target_frame = char(lower(string(params.target_frame)));
problem_data.vinf_target_rth = params.vinf_target_rth_km_s(:);
meta.target_frame = problem_data.target_frame;
meta.vinf_target_rth = problem_data.vinf_target_rth;

% Refresh the fixed-target vector used by optimization mode.
if strcmp(problem_data.mode, 'opt')
    problem_data.vinf_target = vInfPlusEarth(:);
    meta.vinf_target = vInfPlusEarth(:);
end

% Refresh the magnitude seed and target direction used by reachability mode.
if strcmp(problem_data.mode, 'reach')
    vinf_mag_seed = norm(vInfPlusEarth);
    if strcmp(problem_data.target_frame, 'earth_sun_rth')
        vinf_mag_seed = norm(problem_data.vinf_target_rth);
    else
        problem_data.vinf_dir = vInfPlusEarth(:) / vinf_mag_seed;
        meta.vinf_dir = problem_data.vinf_dir;
    end
    if isfield(problem_data, 'vinf_mag_idx') && ~isempty(problem_data.vinf_mag_idx)
        y0(problem_data.vinf_mag_idx) = vinf_mag_seed;
        meta.vInfMag0 = vinf_mag_seed;
    end
    if isfield(opts, 'dv_budget') && ~isempty(opts.dv_budget)
        problem_data.dv_budget = opts.dv_budget;
        meta.dv_budget = opts.dv_budget;
    end
end
end

function [y0, problem_data, meta, lb, ub] = local_update_continuation_time_lock(y0, problem_data, meta, lb, ub, initial_geometry_time)
%LOCAL_UPDATE_CONTINUATION_TIME_LOCK Update the enforced departure-time lock without rebuilding the NLP structure.

% Keep the continuation semantics narrow: only the prescribed initial-time lock is refreshed.
problem_data.initial_geometry_time = initial_geometry_time;
problem_data.t_node1_fixed = initial_geometry_time;
meta.initial_geometry_time = initial_geometry_time;
meta.t_node1_fixed = initial_geometry_time;

if strcmp(problem_data.time_variable_mode, 'propagation_times')
    problem_data.t0_anchor = initial_geometry_time;
    meta.t0_anchor = initial_geometry_time;
    if ~isempty(problem_data.t1_idx)
        y0(problem_data.t1_idx) = initial_geometry_time;
        lb(problem_data.t1_idx) = initial_geometry_time;
        ub(problem_data.t1_idx) = initial_geometry_time;
    end
else
    problem_data.t0_anchor = initial_geometry_time;
    meta.t0_anchor = initial_geometry_time;
    if ~isempty(problem_data.time_idx)
        y0(problem_data.time_idx(1)) = initial_geometry_time;
        lb(problem_data.time_idx(1)) = initial_geometry_time;
        ub(problem_data.time_idx(1)) = initial_geometry_time;
    end
end
end

function seed_time_info = local_build_seed_time_info_from_result(base_result)
%LOCAL_BUILD_SEED_TIME_INFO_FROM_RESULT Rebuild the compact seed-time package from one previous solver result.

seed_time_info = struct();
seed_time_info.mode = base_result.time_mode;
seed_time_info.seed_tof_lower_bound_days = base_result.seed_tof_lower_bound_days;
seed_time_info.seed_tof_exact_days = base_result.seed_tof_exact_days;
seed_time_info.seed_tof_cap_days = base_result.seed_tof_cap_days;
seed_time_info.seed_tof_keep = base_result.seed_tof_keep;
seed_time_info.seed_flyby_inbound_days = base_result.seed_flyby_inbound_days;
seed_time_info.max_tof_relax_factor = base_result.max_tof_relax_factor;
seed_time_info.active_max_tof_days = base_result.active_max_tof_days;
seed_time_info.max_tof_relaxed = base_result.max_tof_relaxed;
seed_time_info.summary_table = base_result.seed_time_summary;
seed_time_info.altitude_scan = base_result.seed_time_altitude_scan;
end

function traj = local_anchor_operational_trajectory_time(traj, params, phase_policy)
%LOCAL_ANCHOR_OPERATIONAL_TRAJECTORY_TIME Anchor operational seed times to the physical epoch constraint before node reconstruction.

% Standard RTH seeds carry a candidate Earth-perigee epoch from the patcher.
if ~isfield(params, 'operational_mode') || strcmpi(string(params.operational_mode), "standard")
    if strcmpi(string(params.target_frame), "earth_sun_rth") && isfield(traj, 'seed_earth_epoch_nd')
        traj = local_shift_trajectory_times(traj, traj.seed_earth_epoch_nd - traj.tTot(traj.earth_peri_idx));
    end
    return
end

% Initial-state mode must start from the epoch implied by the onboard phase.
if strcmpi(string(params.operational_mode), "initial_state")
    departure_epoch = params.initial_epoch_time + phase_policy.phase_seed_nlp - params.initial_rpo_phase;
    time_shift = departure_epoch - traj.tTot(1);
else
    % Fixed-arrival mode shifts the reconstructed seed so Earth perigee lands at the required epoch.
    time_shift = params.final_epoch_time - traj.tTot(traj.earth_peri_idx);
end

% Apply one absolute epoch shift while preserving the reconstructed transfer geometry.
traj = local_shift_trajectory_times(traj, time_shift);
end

function traj = local_shift_trajectory_times(traj, time_shift)
%LOCAL_SHIFT_TRAJECTORY_TIMES Shift every stored trajectory timeline by one absolute epoch offset.

traj.t0 = traj.t0 + time_shift;
traj.tf = traj.tf + time_shift;
traj.tt = traj.tt + time_shift;
traj.ttFB_true = traj.ttFB_true + time_shift;
traj.tSOI = traj.tSOI + time_shift;
traj.t = traj.t + time_shift;
traj.tTot = traj.tTot + time_shift;
traj.tB = traj.tB + time_shift;
traj.tF = traj.tF + time_shift;
traj.ttFB = traj.ttFB + time_shift;
end

function x1_post = local_build_phase_consistent_first_node(y0, phase_theta_idx, x0_rpo_base, I_case2b, mu)
%LOCAL_BUILD_PHASE_CONSISTENT_FIRST_NODE Build the first post-impulse node from theta_rpo and the node-1 impulse seed.

% Propagate the base RPO seed to the active departure phase.
theta_rpo = y0(phase_theta_idx);
x1_post = local_build_phase_consistent_first_node_from_theta(theta_rpo, y0, x0_rpo_base, I_case2b, mu);
end

function x1_post = local_build_phase_consistent_first_node_from_theta(theta_rpo, y0, x0_rpo_base, I_case2b, mu)
%LOCAL_BUILD_PHASE_CONSISTENT_FIRST_NODE_FROM_THETA Build the first post-impulse node from one supplied RPO departure phase.

% Propagate the base RPO seed to the selected departure phase.
x_ini = propagate(x0_rpo_base(:), 0, theta_rpo, mu, '3D', 'CR3BP', 'false');

% Decode the node-1 regularized impulse and map it back to delta-v.
u1 = y0(7);
w1 = y0(8);
s1 = y0(9);
is_case2b = ~isempty(intersect(I_case2b(:), 1));
if ~is_case2b
    dv1 = [u1^2 - w1^2 - s1^2; 2 * u1 * w1; 2 * u1 * s1];
else
    dv1 = [2 * u1 * w1; u1^2 - w1^2 - s1^2; 2 * u1 * s1];
end

% Apply the departure impulse to the propagated phase-dependent boundary state.
x1_post = x_ini(:);
x1_post(4:6) = x1_post(4:6) + dv1(:);
end

function grad_results = run_gradient_checks_phase(mode, y0, problem_data, opts)
%RUN_GRADIENT_CHECKS_PHASE Run MATLAB gradient checks on the phase-aware objective and constraints.

% Build the phase-aware constraint handle chosen by the requested backend.
grad_opts = optimoptions("fmincon", FiniteDifferenceType="central");
nonlcon = get_constraint_handle_phase(mode, problem_data, opts.jacobian_backend);

% Check the objective and constraint gradients in the active mode.
if mode == "opt"
    [validObj, errObj] = checkGradients(@(y) objfun_CR3BP_regularized_phase(y, problem_data), y0, grad_opts, Display="on");
    [validCon, errCon] = checkGradients(nonlcon, y0, grad_opts, IsConstraint=true, Display="on");
else
    [validObj, errObj] = checkGradients(@(y) objfun_CR3BP_reachability_phase(y, problem_data), y0, grad_opts, Display="on");
    [validCon, errCon] = checkGradients(nonlcon, y0, grad_opts, IsConstraint=true, Display="on");
end

% Pack the gradient-check summary in the same format as the legacy branch.
grad_results = struct('validObj', validObj, 'errObj', errObj, ...
    'validCon', validCon, 'errCon', errCon);
end

function solve_out = run_solver_dispatch(mode, y0, problem_data, opts, phase_policy, lb, ub, nodes, impulses, constants, meta)
%RUN_SOLVER_DISPATCH Dispatch the setup to the legacy or phase-aware solver branch.

% Route the normalized IPOPT v1 profile through the CasADi backend.
if opts.solver_backend == "casadi_ipopt"
    solve_out = run_solver_casadi_ipopt(mode, y0, problem_data, opts, lb, ub, nodes, impulses, constants, meta);
    return
end

% Keep the original solver path for seed_phase and switch only the new branches.
if phase_policy.use_phase_variable
    solve_out = run_solver_phase(mode, y0, problem_data, opts, lb, ub, nodes, impulses, constants, meta);
else
    solve_out = run_solver(mode, y0, problem_data, opts, lb, ub, nodes, impulses, constants, meta);
end
end

function solve_out = run_solver_casadi_ipopt(mode, y0, problem_data, opts, lb, ub, nodes, impulses, constants, meta)
%RUN_SOLVER_CASADI_IPOPT Run the normalized initial-state CR3BP profile with CasADi/IPOPT.

% Build the symbolic NLP from the existing packed setup and bounds.
local_progress_log(opts.progress_logs, 'RPO2DS', 'building CasADi NLP; nvars %d, nodes %d', numel(y0), nodes.n_nodes);
problem_definition = build_casadi_rpo2ds_initial_state_problem(y0, lb, ub, problem_data, opts);
local_progress_log(opts.progress_logs, 'RPO2DS', 'CasADi NLP built; equalities %d, inequalities %d', problem_definition.num_equalities, problem_definition.num_inequalities);

% Run the configured IPOPT stages and reuse the existing plot hooks.
local_progress_log(opts.progress_logs, 'RPO2DS', 'calling IPOPT two-stage solver');
solve_out = solve_casadi_ipopt_two_stage(problem_definition, y0, opts);
local_progress_log(opts.progress_logs, 'RPO2DS', 'IPOPT two-stage solver returned; objective %.6g', solve_out.J_opt);
solve_out.stage1_plot_results = local_plot_stage1_solution(solve_out.y_stage1, mode, nodes, impulses, constants, meta, opts);
solve_out.casadi_problem = problem_definition.solver_configuration;
end

function grad_results = run_gradient_checks(mode, y0, problem_data, opts)
%RUN_GRADIENT_CHECKS Compare objective/constraint derivatives at the fixed-phase seed.

grad_opts = optimoptions("fmincon", FiniteDifferenceType="central");
nonlcon = get_constraint_handle(mode, problem_data, opts.jacobian_backend);

if mode == "opt"
    [validObj, errObj] = checkGradients(@(y) objfun_CR3BP_regularized(y, problem_data), y0, grad_opts, Display="on");
    [validCon, errCon] = checkGradients(nonlcon, y0, grad_opts, IsConstraint=true, Display="on");
else
    [validObj, errObj] = checkGradients(@(y) objfun_CR3BP_reachability(y, problem_data), y0, grad_opts, Display="on");
    [validCon, errCon] = checkGradients(nonlcon, y0, grad_opts, IsConstraint=true, Display="on");
end

grad_results = struct('validObj', validObj, 'errObj', errObj, ...
    'validCon', validCon, 'errCon', errCon);
end

function solve_out = run_solver(mode, y0, problem_data, opts, lb, ub, nodes, impulses, constants, meta)
%RUN_SOLVER Run the fixed-phase feasibility and objective solves.
%
% PROTOTYPE:
%   solve_out = run_solver(mode, y0, problem_data, opts, lb, ub, ...
%       nodes, impulses, constants, meta)
%
% DESCRIPTION:
%   Select opt/reach objectives and dense/sparse constraints. Optionally solve
%   zero-objective feasibility, then warm-start the objective stage with the
%   configured iteration and tolerance limits.
%
% INPUT:
%   mode / y0 / problem_data  Objective mode, packed seed, and NLP model.
%   opts / lb / ub            Stage options and variable bound columns.
%   nodes / impulses          Node sampling and impulse-branch metadata.
%   constants / meta          DU/TU/VU scales and packed-variable conventions.
%                               See the public solver header for units/layout.
%
% OUTPUT:
%   solve_out                 Decision vectors, objective, exit flags,
%                               stage-specific diagnostics, and solve options.

nonlcon = get_constraint_handle(mode, problem_data, opts.jacobian_backend);
opt_stage1 = optimoptions('fmincon', ...
    'Display', opts.solver_display, ...
    'Algorithm', opts.stage1_algorithm, ...
    'MaxIter',  opts.max_iter_stage1, ...
    'MaxFunEvals', 1e6, ...
    'ConstraintTolerance', opts.constraint_tolerance, ...
    'TolFun', 1e-6, ...
    'OptimalityTolerance', opts.optimality_tolerance, ...
    'StepTolerance', 1e-6, ...
    'GradObj', 'on', ...
    'GradConstr', 'on', ...
    'DerivativeCheck', 'off', ...
    'FinDiffType', 'central');
opt_stage2 = optimoptions('fmincon', ...
    'Display', opts.solver_display, ...
    'Algorithm', opts.stage2_algorithm, ...
    'MaxIter',  opts.max_iter_stage2, ...
    'MaxFunEvals', 1e6, ...
    'ConstraintTolerance', opts.constraint_tolerance, ...
    'TolFun', 1e-6, ...
    'OptimalityTolerance', opts.optimality_tolerance, ...
    'StepTolerance', 1e-6, ...
    'GradObj', 'on', ...
    'GradConstr', 'on', ...
    'DerivativeCheck', 'off', ...
    'FinDiffType', 'central');

solve_out = struct();
solve_out.stage1_plot_results = struct();

if mode == "opt"
    if local_should_run_stage1(opts)
        [y_stage1, ~, exitflag1, output1] = fmincon(...
            @(y) obj_zero(y), ...
            y0, ...
            [], [], [], [], lb, ub, ...
            nonlcon, ...
            opt_stage1);
    else
        y_stage1 = y0;
        exitflag1 = NaN;
        output1 = struct('message', 'Stage 1 skipped by user option.');
    end

    solve_out.stage1_plot_results = local_plot_stage1_solution(y_stage1, mode, nodes, impulses, constants, meta, opts);

    [y_opt, J_opt, exitflag2, output2] = fmincon(...
        @(y) objfun_CR3BP_regularized(y, problem_data), ...
        y_stage1, ...
        [], [], [], [], lb, ub, ...
        nonlcon, ...
        opt_stage2);

    solve_out.y_opt = y_opt;
    solve_out.y_stage1 = y_stage1;
    solve_out.J_opt = J_opt;
    solve_out.exitflag = exitflag2;
    solve_out.output = output2;
    solve_out.output_stage1 = output1;
    solve_out.exitflag_stage1 = exitflag1;
else
    if local_should_run_stage1(opts)
        [y_stage1, ~, exitflag1, output1] = fmincon(...
            @(y) obj_zero(y), ...
            y0, ...
            [], [], [], [], lb, ub, ...
            nonlcon, ...
            opt_stage1);
    else
        y_stage1 = y0;
        exitflag1 = NaN;
        output1 = struct('message', 'Stage 1 skipped by user option.');
    end

    solve_out.stage1_plot_results = local_plot_stage1_solution(y_stage1, mode, nodes, impulses, constants, meta, opts);

    [y_opt, J_opt, exitflag2, output2] = fmincon(...
        @(y) objfun_CR3BP_reachability(y, problem_data), ...
        y_stage1, ...
        [], [], [], [], lb, ub, ...
        nonlcon, ...
        opt_stage2);

    solve_out.y_opt = y_opt;
    solve_out.y_stage1 = y_stage1;
    solve_out.J_opt = J_opt;
    solve_out.exitflag = 0; % preserve legacy behavior
    solve_out.output = {output1, output2};
    solve_out.output_stage1 = output1;
    solve_out.output_stage2 = output2;
    solve_out.exitflag_stage1 = exitflag1;
    solve_out.exitflag_stage2 = exitflag2;
end

solve_out.options_stage1 = opt_stage1;
solve_out.options_stage2 = opt_stage2;
solve_out.options = opt_stage2;
end

function solve_out = run_solver_phase(mode, y0, problem_data, opts, lb, ub, nodes, impulses, constants, meta)
%RUN_SOLVER_PHASE Run the phase-aware feasibility and objective solves.
%
% PROTOTYPE:
%   solve_out = run_solver_phase(mode, y0, problem_data, opts, lb, ub, ...
%       nodes, impulses, constants, meta)
%
% DESCRIPTION:
%   Select opt/reach objectives and dense/sparse constraints. Optionally solve
%   zero-objective feasibility, then warm-start the objective stage with the
%   configured iteration and tolerance limits.
%
% INPUT:
%   mode / y0 / problem_data  Objective mode, packed seed, and NLP model.
%   opts / lb / ub            Stage options and variable bound columns.
%   nodes / impulses          Node sampling and impulse-branch metadata.
%   constants / meta          DU/TU/VU scales and packed-variable conventions.
%                               See the public solver header for units/layout.
%
% OUTPUT:
%   solve_out                 Decision vectors, objective, exit flags,
%                               stage-specific diagnostics, and solve options.

nonlcon = get_constraint_handle_phase(mode, problem_data, opts.jacobian_backend);
opt_stage1 = optimoptions('fmincon', ...
    'Display', opts.solver_display, ...
    'Algorithm', opts.stage1_algorithm, ...
    'MaxIter', opts.max_iter_stage1, ...
    'MaxFunEvals', 1e6, ...
    'ConstraintTolerance', opts.constraint_tolerance, ...
    'TolFun', 1e-6, ...
    'OptimalityTolerance', opts.optimality_tolerance, ...
    'StepTolerance', 1e-6, ...
    'GradObj', 'on', ...
    'GradConstr', 'on', ...
    'DerivativeCheck', 'off', ...
    'FinDiffType', 'central');
opt_stage2 = optimoptions('fmincon', ...
    'Display', opts.solver_display, ...
    'Algorithm', opts.stage2_algorithm, ...
    'MaxIter', opts.max_iter_stage2, ...
    'MaxFunEvals', 1e6, ...
    'ConstraintTolerance', opts.constraint_tolerance, ...
    'TolFun', 1e-6, ...
    'OptimalityTolerance', opts.optimality_tolerance, ...
    'StepTolerance', 1e-6, ...
    'GradObj', 'on', ...
    'GradConstr', 'on', ...
    'DerivativeCheck', 'off', ...
    'FinDiffType', 'central');

% Run the same feasibility stage used by the legacy branch.
solve_out = struct();
solve_out.stage1_plot_results = struct();
if mode == "opt"
    if local_should_run_stage1(opts)
        [y_stage1, ~, exitflag1, output1] = fmincon( ...
            @(y) obj_zero(y), ...
            y0, [], [], [], [], lb, ub, ...
            nonlcon, ...
            opt_stage1);
    else
        y_stage1 = y0;
        exitflag1 = NaN;
        output1 = struct('message', 'Stage 1 skipped by user option.');
    end

    solve_out.stage1_plot_results = local_plot_stage1_solution(y_stage1, mode, nodes, impulses, constants, meta, opts);

    [y_opt, J_opt, exitflag2, output2] = fmincon( ...
        @(y) objfun_CR3BP_regularized_phase(y, problem_data), ...
        y_stage1, [], [], [], [], lb, ub, ...
        nonlcon, ...
        opt_stage2);

    solve_out.y_opt = y_opt;
    solve_out.y_stage1 = y_stage1;
    solve_out.J_opt = J_opt;
    solve_out.exitflag = exitflag2;
    solve_out.output = output2;
    solve_out.output_stage1 = output1;
    solve_out.output_stage2 = output2;
    solve_out.exitflag_stage1 = exitflag1;
    solve_out.exitflag_stage2 = exitflag2;
else
    if local_should_run_stage1(opts)
        [y_stage1, ~, exitflag1, output1] = fmincon( ...
            @(y) obj_zero(y), ...
            y0, [], [], [], [], lb, ub, ...
            nonlcon, ...
            opt_stage1);
    else
        y_stage1 = y0;
        exitflag1 = NaN;
        output1 = struct('message', 'Stage 1 skipped by user option.');
    end

    solve_out.stage1_plot_results = local_plot_stage1_solution(y_stage1, mode, nodes, impulses, constants, meta, opts);

    [y_opt, J_opt, exitflag2, output2] = fmincon( ...
        @(y) objfun_CR3BP_reachability_phase(y, problem_data), ...
        y_stage1, [], [], [], [], lb, ub, ...
        nonlcon, ...
        opt_stage2);

    solve_out.y_opt = y_opt;
    solve_out.y_stage1 = y_stage1;
    solve_out.J_opt = J_opt;
    solve_out.exitflag = 0;
    solve_out.output = {output1, output2};
    solve_out.output_stage1 = output1;
    solve_out.output_stage2 = output2;
    solve_out.exitflag_stage1 = exitflag1;
    solve_out.exitflag_stage2 = exitflag2;
end

solve_out.options_stage1 = opt_stage1;
solve_out.options_stage2 = opt_stage2;
solve_out.options = opt_stage2;
end

function y0_out = local_apply_y0_override(y0_in, opts)
%LOCAL_APPLY_Y0_OVERRIDE Replace the internally built decision-vector seed when the caller provides one valid warm start.

y0_out = y0_in;
if ~isfield(opts, 'y0_override') || isempty(opts.y0_override)
    return
end

y0_override = opts.y0_override(:);
if numel(y0_override) ~= numel(y0_in)
    error('RPO2DS_transfer:InvalidY0OverrideLength', ...
        'opts.y0_override must have the same length as the internally generated decision vector.');
end

y0_out = y0_override;
end

function run_stage1 = local_should_run_stage1(opts)
%LOCAL_SHOULD_RUN_STAGE1 Read the stage-1 toggle and default to true when the caller does not override it.

run_stage1 = true;
if isfield(opts, 'run_stage1') && ~isempty(opts.run_stage1)
    run_stage1 = logical(opts.run_stage1);
end
end

function stage1_plot_results = local_plot_stage1_solution(y_stage1, mode, nodes, impulses, constants, meta, opts)
%LOCAL_PLOT_STAGE1_SOLUTION Plot the stage-1 solution immediately after the first fmincon solve when requested.

stage1_plot_results = struct();
if ~local_should_run_stage1(opts)
    return
end
if isempty(y_stage1)
    return
end
if ~opts.graph.objzero_rot && ~opts.graph.objzero_inertial
    return
end

stage1_post = build_preopt_post(mode, y_stage1, nodes, impulses, constants, meta, opts);
if opts.graph.objzero_rot
    stage1_plot_results.objzero_rot = plot_objzero_rot(stage1_post, constants);
end
if opts.graph.objzero_inertial
    stage1_plot_results.objzero_inertial = plot_objzero_inertial(stage1_post, constants);
end
drawnow;
end

function nonlcon = get_constraint_handle_phase(mode, problem_data, jacobian_backend)
%GET_CONSTRAINT_HANDLE_PHASE Build the phase-aware constraint handle in dense or sparse form.

% Select the requested dense or sparse wrapper in the active solve mode.
if mode == "opt"
    if jacobian_backend == "sparse"
        nonlcon = @(y) constr_CR3BP_regularized_phase_sparse(y, problem_data);
    else
        nonlcon = @(y) constr_CR3BP_regularized_phase(y, problem_data);
    end
else
    if jacobian_backend == "sparse"
        nonlcon = @(y) constr_CR3BP_reachability_phase_sparse(y, problem_data);
    else
        nonlcon = @(y) constr_CR3BP_reachability_phase(y, problem_data);
    end
end
end

function nonlcon = get_constraint_handle(mode, problem_data, jacobian_backend)
%GET_CONSTRAINT_HANDLE Select opt/reach constraints in the requested dense or sparse storage.

if mode == "opt"
    if jacobian_backend == "sparse"
        nonlcon = @(y) constr_CR3BP_regularized_sparse(y, problem_data);
    else
        nonlcon = @(y) constr_CR3BP_regularized(y, problem_data);
    end
else
    if jacobian_backend == "sparse"
        nonlcon = @(y) constr_CR3BP_reachability_sparse(y, problem_data);
    else
        nonlcon = @(y) constr_CR3BP_reachability(y, problem_data);
    end
end
end

function post = local_build_preopt_setup(mode, y_seed, nodes, impulses, constants, meta, opts)
%LOCAL_BUILD_PREOPT_SETUP Build the setup diagnostic payload without propagating unless plots need it.

if local_needs_full_preopt_post(opts)
    post = build_preopt_post(mode, y_seed, nodes, impulses, constants, meta, opts);
else
    post = local_build_light_preopt_post(mode, y_seed, nodes, impulses, constants, meta);
end
end

function needs_full = local_needs_full_preopt_post(opts)
%LOCAL_NEEDS_FULL_PREOPT_POST Check whether setup diagnostics need full propagated trajectories.

needs_full = (~isfield(opts, 'setup_only') || ~opts.setup_only) || ...
    (isfield(opts, 'graph') && ...
    (opts.graph.initial_guess_fmincon_rot || opts.graph.initial_guess_fmincon_inertial));
end

function post = local_ensure_full_preopt_post(post, mode, y_seed, nodes, impulses, constants, meta, opts)
%LOCAL_ENSURE_FULL_PREOPT_POST Upgrade a lightweight setup diagnostic payload when a trajectory plot needs it.

if isfield(post, 'is_lightweight') && post.is_lightweight
    post = build_preopt_post(mode, y_seed, nodes, impulses, constants, meta, opts);
end
end

function post = local_build_light_preopt_post(mode, y_seed, nodes, impulses, constants, meta)
%LOCAL_BUILD_LIGHT_PREOPT_POST Rebuild node-level setup diagnostics without segment propagation.

if isfield(meta, 'use_phase_variable') && meta.use_phase_variable
    post = local_build_light_preopt_post_phase(mode, y_seed, nodes, impulses, constants, meta);
    return
end

post = local_build_light_preopt_post_core(mode, y_seed, nodes, impulses, constants, meta);
end

function post = local_build_light_preopt_post_phase(mode, y_seed, nodes, impulses, constants, meta)
%LOCAL_BUILD_LIGHT_PREOPT_POST_PHASE Rebuild lightweight setup diagnostics for the phase-aware decision vector.

if isfield(meta, 'has_phase_variable') && ~meta.has_phase_variable
    y_legacy = y_seed;
    t_nodes_phase = local_decode_time_nodes(y_seed, meta);
    phase_theta = meta.initial_rpo_phase + t_nodes_phase(1) - meta.initial_epoch_time;
else
    phase_theta = y_seed(meta.phase_theta_idx);
    if mode == "opt"
        y_legacy = y_seed(1:end-1);
    else
        y_legacy = [y_seed(1:meta.phase_theta_idx-1); y_seed(meta.phase_theta_idx+1:end)];
    end
end

post = local_build_light_preopt_post_core(mode, y_legacy, nodes, impulses, constants, meta);
post.phase_rpo = phase_theta;
post.phase_seed = meta.phase_seed;
post.phase_request = meta.phase_request;
phase_theta_geometry = phase_theta;
if isfield(meta, 'phase_origin') && isfield(meta, 'phase_period') && ~isempty(meta.phase_origin)
    phase_theta_geometry = local_wrap_phase_value(phase_theta, meta.phase_origin, meta.phase_period);
end
post.phase_rpo_wrapped = phase_theta_geometry;
post.x0_rpo_phase = [];
end

function post = local_build_light_preopt_post_core(mode, y_seed, nodes, impulses, constants, meta)
%LOCAL_BUILD_LIGHT_PREOPT_POST_CORE Decode node-level setup diagnostics without integrating between nodes.

mu_EM = constants.mu_EM;
VU = constants.VU;
n_nodes = nodes.n_nodes;
t_nodes = local_decode_time_nodes(y_seed, meta);
X_nodes = reshape(y_seed(1:9*n_nodes), [9, n_nodes]);
xx_nodes = X_nodes(1:6,:)';
U_nodes = X_nodes(7:9,:)';
vInfMag = [];
if mode == "reach"
    vInfMag = y_seed(end);
end

xx_nodes_in = rot2in(xx_nodes, t_nodes, mu_EM);
xx_nodes_in(:,1:3) = xx_nodes_in(:,1:3) * constants.DU;
xx_nodes_in(:,4:6) = xx_nodes_in(:,4:6) * constants.VU;
[DV_mag_opt, DVvec_opt_rot, DVvec_opt_in] = compute_dv_vectors(U_nodes, impulses.I_case2b, xx_nodes, t_nodes, mu_EM, VU);

post = struct();
post.is_lightweight = true;
post.dv_opt = (1000 * VU) * sum(DV_mag_opt);
post.tof_opt = (constants.TU / 86400) * (t_nodes(end) - t_nodes(1));
post.t_nodes = t_nodes;
post.X_nodes = X_nodes;
post.xx_nodes = xx_nodes;
post.U_nodes = U_nodes;
post.xx_traj_rot = [];
post.t_traj = [];
post.STM_vec = cell(n_nodes - 1, 1);
post.xx_traj_in = [];
post.xx_nodes_in = xx_nodes_in;
post.DV_mag_opt = DV_mag_opt;
post.DVvec_opt_rot = DVvec_opt_rot;
post.DVvec_opt_in = DVvec_opt_in;
post.vInfMag = vInfMag;
end

function post = build_preopt_post(mode, y_seed, nodes, impulses, constants, meta, opts)
%BUILD_PREOPT_POST Rebuild one NLP decision vector into the shared trajectory struct.

% Silence compact logs while reconstructing auxiliary plot trajectories.
opts_quiet = opts;
opts_quiet.print_logs = false;

% Reuse the phase-aware or fixed-phase post-processing branch as needed.
if isfield(meta, 'use_phase_variable') && meta.use_phase_variable
    post = postprocess_solution_phase(mode, y_seed, nodes, impulses, constants, meta, opts_quiet);
else
    post = postprocess_solution(mode, y_seed, nodes, impulses, constants, meta, opts_quiet);
end
post.is_lightweight = false;
end

function post = postprocess_solution(mode, y_opt, nodes, impulses, constants, meta, opts)
%POSTPROCESS_SOLUTION Decode nodes, propagate segments, and recover burns.
%
% PROTOTYPE:
%   post = postprocess_solution(mode, y_opt, nodes, impulses, constants, meta, opts)
%
% DESCRIPTION:
%   Reconstruct absolute epochs, map regularized U to velocity jumps,
%   propagate segment histories, and transform them to Earth-centered
%   inertial axes.
%
% INPUT:
%   mode / y_opt              Active opt/reach mode and optimized vector.
%   nodes / impulses / meta   Shooting layout, impulse branches, and decoding.
%   constants / opts          Physical scales and requested diagnostic plots.
%
% OUTPUT:
%   post                      Node/segment histories, impulses, physical
%                               delta-v [m/s] and TOF [days]. Rotating
%                               coordinates use DU/VU; inertial positions km.

mu_EM = constants.mu_EM;
DU = constants.DU;
VU = constants.VU;

n_nodes = nodes.n_nodes;

% Reconstruct the absolute node times from the active time-variable layout.
t_nodes = local_decode_time_nodes(y_opt, meta);

% Extract the optimized state and impulse variables
X_nodes = reshape(y_opt(1:9*n_nodes), [9, n_nodes]);
xx_nodes = X_nodes(1:6,:)';
U_nodes = X_nodes(7:9,:)';

% Extract magnitude of vInf for reachability
vInfMag = [];
if mode == "reach"
    vInfMag = y_opt(end);
end

% Compute trajectory in rotating frame
xx_traj_rot = [];
t_traj = [];
STM_vec = cell(n_nodes-1,1);

dimensions = '3D';
framework = 'CR3BP';
stm_computation = 'true';

for j = 1:n_nodes-1
    [~,~,STM_hist, x_temp, t_temp] = propagate(xx_nodes(j,:), t_nodes(j), t_nodes(j+1), mu_EM, ...
        dimensions, framework, stm_computation);
    x_temp = x_temp(:,1:6);
    xx_traj_rot = [xx_traj_rot; x_temp];
    t_traj = [t_traj; t_temp];
    STM_vec{j} = STM_hist;
end

% Transform into ECI
xx_traj_in = rot2in(xx_traj_rot, t_traj, mu_EM);
xx_traj_in(:,1:3) = xx_traj_in(:,1:3) * DU;
xx_traj_in(:,4:6) = xx_traj_in(:,4:6) * VU;

xx_nodes_in = rot2in(xx_nodes, t_nodes, mu_EM);
xx_nodes_in(:,1:3) = xx_nodes_in(:,1:3) * DU;
xx_nodes_in(:,4:6) = xx_nodes_in(:,4:6) * VU;

% Compute impulses
[DV_mag_opt, DVvec_opt_rot, DVvec_opt_in] = compute_dv_vectors(U_nodes, impulses.I_case2b, xx_nodes, t_nodes, mu_EM, VU);

% total DV
post.dv_opt = (1000*VU) * sum(DV_mag_opt); % [m/s]

% total TOF
post.tof_opt = (constants.TU/86400) * (t_nodes(end)-t_nodes(1)); % [days]

% Verify distance from Moon and Earth
[~, ~, is_below_minM] = distance_from_point(xx_traj_rot, 1-mu_EM, constants.rminMoon_DU);
[~, ~, is_below_minE] = distance_from_point(xx_traj_rot, -mu_EM, constants.rminEarth_DU);

% Display results
if opts.print_logs
    fprintf('Total Delta V: %.4f m/s\n', post.dv_opt);
    fprintf('Total Time of Flight: %.4f days\n', post.tof_opt);
    if mode == "reach"
        fprintf('Reachable vInf magnitude: %.6f km/s\n', vInfMag);
    end
    if is_below_minM == 0
        fprintf('Moon distance limit respected\n');
    else
        fprintf('Moon distance limit not respected\n');
    end
    if is_below_minE == 0
        fprintf('Earth distance limit respected\n');
    else
        fprintf('Earth distance limit not respected\n');
    end
end

post.t_nodes = t_nodes;
post.X_nodes = X_nodes;
post.xx_nodes = xx_nodes;
post.U_nodes = U_nodes;
post.xx_traj_rot = xx_traj_rot;
post.t_traj = t_traj;
post.STM_vec = STM_vec;
post.xx_traj_in = xx_traj_in;
post.xx_nodes_in = xx_nodes_in;
post.DV_mag_opt = DV_mag_opt;
post.DVvec_opt_rot = DVvec_opt_rot;
post.DVvec_opt_in = DVvec_opt_in;
post.vInfMag = vInfMag;

end

function post = postprocess_solution_phase(mode, y_opt, nodes, impulses, constants, meta, opts)
%POSTPROCESS_SOLUTION_PHASE Decode phase-aware nodes, propagate segments, and recover burns.
%
% PROTOTYPE:
%   post = postprocess_solution_phase(mode, y_opt, nodes, impulses, constants, meta, opts)
%
% DESCRIPTION:
%   Recover the active RPO departure phase/state and time layout, then reuse
%   the shared trajectory/impulse reconstruction.
%
% INPUT:
%   mode / y_opt              Active opt/reach mode and optimized vector.
%   nodes / impulses / meta   Shooting layout, impulse branches, and decoding.
%   constants / opts          Physical scales and requested diagnostic plots.
%
% OUTPUT:
%   post                      Node/segment histories, impulses, physical
%                               delta-v [m/s] and TOF [days]. Rotating
%                               coordinates use DU/VU; inertial positions km.

if isfield(meta, 'has_phase_variable') && ~meta.has_phase_variable
    y_legacy = y_opt;
    t_nodes_phase = local_decode_time_nodes(y_opt, meta);
    phase_theta = meta.initial_rpo_phase + t_nodes_phase(1) - meta.initial_epoch_time;
else
    phase_theta = y_opt(meta.phase_theta_idx);
    if mode == "opt"
        y_legacy = y_opt(1:end-1);
    else
        y_legacy = [y_opt(1:meta.phase_theta_idx-1); y_opt(meta.phase_theta_idx+1:end)];
    end
end

% Reuse the legacy trajectory and delta-v extraction on the state/time variables.
post = postprocess_solution(mode, y_legacy, nodes, impulses, constants, meta, opts);

% Store the phase diagnostics required by the new architecture.
post.phase_rpo = phase_theta;
post.phase_seed = meta.phase_seed;
post.phase_request = meta.phase_request;
phase_theta_geometry = phase_theta;
if isfield(meta, 'phase_origin') && isfield(meta, 'phase_period') && ~isempty(meta.phase_origin)
    phase_theta_geometry = local_wrap_phase_value(phase_theta, meta.phase_origin, meta.phase_period);
end
post.phase_rpo_wrapped = phase_theta_geometry;
post.x0_rpo_phase = propagate(constants.xx0RPOBase, 0, phase_theta_geometry, constants.mu_EM, '3D', 'CR3BP', 'false');

% Print the optimized phase only when compact logging is enabled.
if opts.print_logs
    fprintf('RPO departure phase: %.6f\n', post.phase_rpo);
end
end

function t_nodes = local_decode_time_nodes(y, meta)
%LOCAL_DECODE_TIME_NODES Rebuild the absolute node times from the packed time variables.

if strcmp(meta.time_variable_mode, 'propagation_times')
    dt_nodes = y(meta.dt_idx(:)); % Read the optimized segment durations
    if isfield(meta, 'time_anchor_mode') && strcmp(meta.time_anchor_mode, 'final')
        t_nodes = meta.t0_anchor - [flipud(cumsum(flipud(dt_nodes))); 0]; % Rebuild times backward from fixed arrival
    elseif ~isempty(meta.t1_idx)
        t1 = y(meta.t1_idx); % Read the optimized departure time when it is part of the NLP
        t_nodes = t1 + [0; cumsum(dt_nodes)]; % Rebuild the absolute node times from the active departure time
    else
        t1 = meta.t0_anchor; % Use the fixed departure time anchor when geometry locks t(1)
        t_nodes = t1 + [0; cumsum(dt_nodes)]; % Rebuild the absolute node times from the active departure time
    end
else
    t_nodes = y(meta.time_idx(:)); % Read the legacy absolute node times directly
end
end

function phase_wrapped = local_wrap_phase_value(phase_raw, phase_origin, phase_period)
%LOCAL_WRAP_PHASE_VALUE Wrap one RPO phase onto the geometry interval used by the database.

% Keep non-periodic values unchanged when the period metadata is absent.
phase_wrapped = phase_raw;
if isempty(phase_period) || ~isfinite(phase_period) || phase_period <= 0
    return
end

% Shift by whole periods while preserving the database phase origin.
phase_wrapped = phase_origin + mod(phase_raw - phase_origin, phase_period);
end

function [DV_mag_opt, DVvec_opt_rot, DVvec_opt_in] = compute_dv_vectors(U_nodes, I_case2b, xx_nodes, t_nodes, mu_EM, VU)
% DV magnitude
DV_mag_opt = U_nodes(:,1).^2 + U_nodes(:,2).^2 + U_nodes(:,3).^2;

% DV vector (rotating frame)
DVx_opt = U_nodes(:,1).^2 - U_nodes(:,2).^2 - U_nodes(:,3).^2;
DVy_opt = 2 * U_nodes(:,1) .* U_nodes(:,2);
DVz_opt = 2 * U_nodes(:,1) .* U_nodes(:,3);

if ~isempty(I_case2b)
    DVx_opt(I_case2b) = 2 * U_nodes(I_case2b,1) .* U_nodes(I_case2b,2);
    DVy_opt(I_case2b) = U_nodes(I_case2b,1).^2 - U_nodes(I_case2b,2).^2 - U_nodes(I_case2b,3).^2;
    DVz_opt(I_case2b) = 2 * U_nodes(I_case2b,1) .* U_nodes(I_case2b,3);
end

DVvec_opt_rot = [DVx_opt, DVy_opt, DVz_opt];

% DV vector inertial
inertial_no_dv = rot2in(xx_nodes, t_nodes, mu_EM);
inertial_with_dv = rot2in([xx_nodes(:,1:3), xx_nodes(:,4:6) + DVvec_opt_rot], t_nodes, mu_EM);
diff = inertial_with_dv - inertial_no_dv;
DVvec_opt_in = diff(:,4:6) * VU * 1000; % [m/s]

end

function plot_ok = plot_flyby(traj, constants)
%PLOT_FLYBY Display the reconstructed lunar flyby in the rotating frame.

figure;
hold on;
grid on;
plot3(traj.xxF(1,1), traj.xxF(1,2), traj.xxF(1,3), 'r.', 'MarkerSize', 20, 'DisplayName','Perilune');
plot3(1-constants.mu_EM, 0, 0, 'k.', 'MarkerSize', 20, 'DisplayName','Moon');
plot3(traj.xxF(:,1), traj.xxF(:,2), traj.xxF(:,3), 'r', 'LineWidth', 2, 'DisplayName','Forward');
plot3(traj.xxB(:,1), traj.xxB(:,2), traj.xxB(:,3), 'b', 'LineWidth', 2, 'DisplayName','Backward');
xlabel('X (DU)', 'Interpreter', 'latex');
ylabel('Y (DU)', 'Interpreter', 'latex');
zlabel('Z (DU)', 'Interpreter', 'latex');
title('Moon swing-by - Adimensional EMB Centered Rotating Frame', 'Interpreter','latex');
set(gcf,"Theme","light");
legend('show', 'interpreter', 'latex');
axis equal;
xlim([0.8 1.05]);
ylim([min(traj.xxFB(:,2)) max(traj.xxFB(:,2))]);
grid on;
hold off;
set(gca,'FontSize',16,'FontName','times new roman');

plot_ok = true;
end

function plot_ok = plot_initial_guess(traj, constants)
%PLOT_INITIAL_GUESS Display the joined patched initial trajectory and its endpoints.

figure;
hold on;
plot3(traj.xx(:,1),traj.xx(:,2),traj.xx(:,3), 'LineWidth',1.5,  'DisplayName', 'RPO-Moon leg');
plot3(traj.xxRot(:,1), traj.xxRot(:,2),traj.xxRot(:,3), 'LineWidth',1.5,  'DisplayName', 'Moon-Earth leg');
theta = linspace(0, 2*pi, 400);
xc = (1 - constants.mu_EM) + constants.rSOI*cos(theta);
yc = constants.rSOI*sin(theta);
zc = - ones(size(theta));
fill3(xc, yc, zc, [1.0  1.0  0.6], 'FaceAlpha', 1, 'EdgeColor', 'none', 'DisplayName', 'Moon SOI');
plot3(traj.xxFB(:,1), traj.xxFB(:,2), traj.xxFB(:,3), 'LineWidth', 1.5, 'DisplayName','Moon Flyby leg');
plot3(traj.xx(1,1),traj.xx(1,2),traj.xx(1,3), 'o', 'Color','b','MarkerSize', 10, 'LineWidth', 3, 'DisplayName','Initial point');
plot3(traj.xxRot(end,1), traj.xxRot(end,2),traj.xxRot(end,3), 'o', 'Color','r','MarkerSize', 10, 'LineWidth', 3, 'DisplayName','Final point');
axis equal;
xlabel('X (DU)', 'Interpreter', 'latex');
ylabel('Y (DU)', 'Interpreter', 'latex');
zlabel('Z (DU)', 'Interpreter', 'latex');
title('Initial Guess Trajectory - Adimensional EMB Centered Rotating Frame', 'Interpreter','latex');
set(gcf,"Theme","light");
legend('show', 'Interpreter', 'Latex');
xlim([min(traj.xxTot(:,1))*1.05 max(traj.xxTot(:,1))*1.75]);
ylim([min(traj.xxTot(:,2)) max(traj.xxTot(:,2))]);
grid on;
hold off;
set(gca,'FontSize',16,'FontName','times new roman');

plot_ok = true;
end

function plot_ok = plot_nodes(traj, nodes, constants)
%PLOT_NODES Overlay the selected multiple-shooting nodes on the seed trajectory.

figure;
hold on;
box on;
grid on;
plot3(traj.xxTot(:,1), traj.xxTot(:,2), traj.xxTot(:,3), '-', 'color', 'r','DisplayName', 'Trajectory');
plot3(nodes.X_nodes_IG(:,1), nodes.X_nodes_IG(:,2), nodes.X_nodes_IG(:,3), 'o', 'color','k', 'DisplayName', 'Nodes');
plot3(1-constants.mu_EM, 0, 0, '.', 'color', 'k', 'MarkerSize', 10, 'DisplayName', 'Moon');
plot3(-constants.mu_EM, 0, 0, '.', 'color', 'b', 'MarkerSize', 15, 'DisplayName', 'Earth');
axis equal;
title('Initial Guess Trajectory with Nodes - Adimensional EMB Centered Rotating Frame', 'Interpreter','latex');
set(gca, 'FontName','Times','FontSize',16);
xlabel('x (DU)', 'Interpreter', 'Latex');
ylabel('y (DU)', 'Interpreter', 'Latex');
zlabel('z (DU)', 'Interpreter', 'Latex');
legend('show', 'Interpreter', 'latex');
xlim([min(traj.xxTot(:,1))*1.05 max(traj.xxTot(:,1))*1.5]);
ylim([min(traj.xxTot(:,2)) max(traj.xxTot(:,2))]);
set(gcf,"Theme","light");

plot_ok = true;
end

function plot_ok = plot_optimized_rot(post, ~)
%PLOT_OPTIMIZED_ROT Plot the optimized rotating-frame trajectory with the shared style.
plot_ok = plot_post_rot(post, 'Optimized Trajectory - Adimensional EMB Centered Rotating Frame', 'Optimized trajectory');
end

function plot_ok = plot_optimized_inertial(post, constants)
%PLOT_OPTIMIZED_INERTIAL Plot the optimized inertial-frame trajectory with the shared style.
plot_ok = plot_post_inertial(post, constants, 'Optimized Trajectory - Earth Centered Inertial Frame', 'Optimized trajectory');
end

function plot_ok = plot_initial_guess_fmincon_rot(post, ~)
%PLOT_INITIAL_GUESS_FMINCON_ROT Plot the NLP initial guess in the rotating frame with the optimized style.
plot_ok = plot_post_rot(post, 'Initial Guess for Fmincon - Adimensional EMB Centered Rotating Frame', 'Initial guess');
end

function plot_ok = plot_initial_guess_fmincon_inertial(post, constants)
%PLOT_INITIAL_GUESS_FMINCON_INERTIAL Plot the NLP initial guess in the inertial frame with the optimized style.
plot_ok = plot_post_inertial(post, constants, 'Initial Guess for Fmincon - Earth Centered Inertial Frame', 'Initial guess');
end

function plot_ok = plot_objzero_rot(post, ~)
%PLOT_OBJZERO_ROT Plot the stage-1 feasibility solution in the rotating frame with the optimized style.
plot_ok = plot_post_rot(post, 'Stage-1 Objzero Trajectory - Adimensional EMB Centered Rotating Frame', 'Stage-1 objzero trajectory');
end

function plot_ok = plot_objzero_inertial(post, constants)
%PLOT_OBJZERO_INERTIAL Plot the stage-1 feasibility solution in the inertial frame with the optimized style.
plot_ok = plot_post_inertial(post, constants, 'Stage-1 Objzero Trajectory - Earth Centered Inertial Frame', 'Stage-1 objzero trajectory');
end

function plot_ok = plot_post_rot(post, title_text, traj_label)
%PLOT_POST_ROT Draw one rotating-frame trajectory with the shared optimized style.

figure;
hold on;
box on;
grid on;
plot3(post.xx_traj_rot(:,1), post.xx_traj_rot(:,2), post.xx_traj_rot(:,3), '-', 'Color', 'k', 'LineWidth', 1, 'DisplayName', traj_label);
plot3(post.xx_nodes(:,1), post.xx_nodes(:,2), post.xx_nodes(:,3), 'o', 'Color', 'k', 'DisplayName', 'Nodes');
quiver3(post.xx_nodes(:,1), post.xx_nodes(:,2), post.xx_nodes(:,3), ...
    post.DVvec_opt_rot(:,1), post.DVvec_opt_rot(:,2), post.DVvec_opt_rot(:,3), ...
    2, 'Color', 'b', 'AutoScale', 'off', 'DisplayName', '$\Delta$V impulses');

title(title_text, 'Interpreter', 'Latex');
set(gca, 'FontName', 'Times', 'FontSize', 15);
xlabel('x (DU)', 'Interpreter', 'Latex');
ylabel('y (DU)', 'Interpreter', 'Latex');
zlabel('z (DU)', 'Interpreter', 'Latex');
axis equal;
local_apply_equal_xy_limits(gca, post.xx_traj_rot(:,1), post.xx_traj_rot(:,2));
legend('show', 'Interpreter', 'latex');
set(gca, 'FontSize', 15, 'FontName', 'times new roman');
set(gcf, "Theme", "light");

plot_ok = true;
end

function plot_ok = plot_post_inertial(post, constants, title_text, traj_label)
%PLOT_POST_INERTIAL Draw one inertial-frame trajectory with the shared optimized style.

% Build the Moon path used by all inertial trajectory plots.
thetaM = linspace(0, 2*pi, 100)';
moonPath = [cos(thetaM), sin(thetaM), zeros(length(thetaM), 1)] * constants.DU;

figure;
hold on;
box on;
grid on;
plot3(post.xx_traj_in(:,1), post.xx_traj_in(:,2), post.xx_traj_in(:,3), '-', 'Color', 'k', 'LineWidth', 1, 'DisplayName', traj_label);
plot3(post.xx_nodes_in(:,1), post.xx_nodes_in(:,2), post.xx_nodes_in(:,3), 'o', 'Color', 'k', 'DisplayName', 'Nodes');
plot3(moonPath(:,1), moonPath(:,2), moonPath(:,3), '--', 'Color', [0.5, 0.5, 0.5], 'DisplayName', 'Moon Path');
quiver3(post.xx_nodes_in(:,1), post.xx_nodes_in(:,2), post.xx_nodes_in(:,3), ...
    post.DVvec_opt_in(:,1), post.DVvec_opt_in(:,2), post.DVvec_opt_in(:,3), ...
    2, 'Color', 'b', 'AutoScale', 'off', 'DisplayName', '$\Delta$V impulses');

axis equal;
local_apply_equal_xy_limits(gca, [post.xx_traj_in(:,1); moonPath(:,1)], [post.xx_traj_in(:,2); moonPath(:,2)]);
title(title_text, 'Interpreter', 'Latex');
set(gca, 'FontName', 'Times', 'FontSize', 15);
xlabel('x (km)', 'Interpreter', 'Latex');
ylabel('y (km)', 'Interpreter', 'Latex');
zlabel('z (km)', 'Interpreter', 'Latex');
legend('show', 'Interpreter', 'latex');
set(gca, 'FontSize', 15, 'FontName', 'times new roman');
set(gcf, "Theme", "light");

plot_ok = true;
end

function local_apply_equal_xy_limits(ax_handle, x_values, y_values)
%LOCAL_APPLY_EQUAL_XY_LIMITS Apply one square XY viewport so rotating and inertial plots share the same visual aspect.

% Use finite plotted coordinates only when computing the visible data bounds.
x_values = x_values(isfinite(x_values));
y_values = y_values(isfinite(y_values));

% Compute a common span so one coordinate direction cannot stretch the figure.
x_center = 0.5 * (min(x_values) + max(x_values));
y_center = 0.5 * (min(y_values) + max(y_values));
x_span = max(x_values) - min(x_values);
y_span = max(y_values) - min(y_values);
xy_span = 1.10 * max([x_span, y_span, eps]);

% Set matching X and Y windows around the plotted data center.
xlim(ax_handle, x_center + 0.5 * xy_span * [-1, 1]);
ylim(ax_handle, y_center + 0.5 * xy_span * [-1, 1]);
pbaspect(ax_handle, [1, 1, 1]);
end

function result = build_result_struct(mode, FullSolution, vInfPlusEarth, xx0_RPO, params, opts, constants, traj, nodes, impulses, meta, solve_out, post, preopt_results, plot_results, exitflag, output, h, guess_diagnostics, seed_time_info, grad_results, setup_data)
%BUILD_RESULT_STRUCT Collect the complete transfer setup, solve, and trajectory result.
%
% PROTOTYPE:
%   result = build_result_struct(mode, FullSolution, vInfPlusEarth, xx0_RPO, ...
%       params, opts, constants, traj, nodes, impulses, meta, solve_out, post, ...
%       preopt_results, plot_results, exitflag, output, h, guess_diagnostics, ...
%       seed_time_info, grad_results, setup_data)
%
% DESCRIPTION:
%   Preserve selected patch data, normalized policies, stage diagnostics,
%   gradient-check results, and plot status with physical output metrics.
%   A returned solver status is not independent numerical acceptance.
%
% INPUT:
%   mode / vInfPlusEarth / xx0_RPO  Solve mode, physical target [km/s],
%                                  and planar rotating seed [DU/VU].
%   FullSolution / h              Patch selection and Earth altitude [km].
%   params / opts / constants     Normalized policies, solver options, scales.
%   traj / nodes / impulses / meta Seed trajectory and NLP decoding payloads.
%   solve_out / exitflag / output  Final solve payload and status.
%   post / preopt_results         Optimized and seed/stage-one geometry.
%   plot_results / grad_results   Optional plotting/derivative diagnostics.
%   guess_diagnostics / seed_time_info / setup_data
%                                 Patch provenance, time caps, reusable NLP setup.
%
% OUTPUT:
%   result  Full solver result. Rotating states use DU/VU and TU; physical
%             delta-v/TOF use m/s/days. Field conventions are described in
%             the public RPO2DS_transfer header.

result = struct();

% Primary outputs
if nargin >= 12 && ~isempty(post)
    result.dv_opt = post.dv_opt;
    result.tof_opt = post.tof_opt;
else
    result.dv_opt = [];
    result.tof_opt = [];
end

if nargin >= 10 && ~isempty(solve_out)
    result.y_opt = solve_out.y_opt;
    result.J_opt = solve_out.J_opt;
    result.exitflag = exitflag;
    result.output = output;
else
    result.y_opt = [];
    result.J_opt = [];
    result.exitflag = [];
    result.output = [];
end

if mode == "reach"
    if nargin >= 12 && ~isempty(post)
        result.vInfMag = post.vInfMag;
    else
        result.vInfMag = [];
    end
end

result.h = h;

% FullSolution
result.FullSolution = FullSolution;

% Rich structured results
result.inputs = struct('mode', mode, 'vInfPlusEarth', vInfPlusEarth, 'xx0_RPO', xx0_RPO, ...
    'params', params, 'opts', opts);

result.constants = constants;
result.initial_guess = traj;
result.nodes = nodes;
result.impulses = impulses;
result.meta = meta;
result.setup = setup_data;
result.initial_guess_search = guess_diagnostics;
result.phase_mode = meta.phase_mode;
result.phase_seed = meta.phase_seed;
result.phase_request = meta.phase_request;
result.phase_shortlist_rows = meta.phase_shortlist_rows;
result.phase_shortlist_phases = meta.phase_shortlist_phases;
result.phase_shortlist_distances = meta.phase_shortlist_distances;
if ~isfield(meta, 'target_frame')
    meta.target_frame = 'eci';
    meta.theta_sun0_mode = 'specified';
    meta.theta_sun0_value = 0;
    meta.theta_sun0_seed = 0;
    meta.theta_sun0_idx = [];
    meta.ws_sun_em_nd = -9.25195985e-1;
    meta.vinf_target_rth = [0; 3; 0];
end
result.target_frame = meta.target_frame;
result.theta_sun0_mode = meta.theta_sun0_mode;
result.theta_sun0_value = meta.theta_sun0_value;
result.theta_sun0_seed = meta.theta_sun0_seed;
result.theta_sun0_idx = meta.theta_sun0_idx;
result.ws_sun_em_nd = meta.ws_sun_em_nd;
result.vinf_target_rth = meta.vinf_target_rth;
result.operational_mode = meta.operational_mode;
result.initial_epoch_time = meta.initial_epoch_time;
result.initial_rpo_phase = meta.initial_rpo_phase;
result.final_epoch_time = meta.final_epoch_time;
result.available_time_nd = meta.available_time_nd;
result.seed_candidate_index = meta.seed_candidate_index;
result.operational_screening_summary = meta.operational_screening_summary;
result.time_mode = meta.time_mode;
result.time_variable_mode = meta.time_variable_mode;
result.geometry_mode = meta.geometry_mode;
result.initial_geometry_time = meta.initial_geometry_time;
result.seed_tof_lower_bound_days = meta.seed_tof_lower_bound_days;
result.seed_tof_exact_days = meta.seed_tof_exact_days;
result.seed_tof_cap_days = meta.seed_tof_cap_days;
result.seed_tof_keep = meta.seed_tof_keep;
result.seed_flyby_inbound_days = meta.seed_flyby_inbound_days;
result.max_tof_relax_factor = meta.max_tof_relax_factor;
result.active_max_tof_days = meta.active_max_tof_days;
result.max_tof_relaxed = meta.max_tof_relaxed;
result.seed_time_summary = seed_time_info.summary_table;
result.seed_time_altitude_scan = seed_time_info.altitude_scan;
result.phase_final = [];
if ~isempty(post) && isfield(post, 'phase_rpo')
    result.phase_final = post.phase_rpo;
elseif ~isempty(setup_data) && isfield(meta, 'phase_source') && strcmp(meta.phase_source, 'initial_time') && isfield(setup_data, 'y0')
    setup_t_nodes_phase = local_decode_time_nodes(setup_data.y0, meta);
    result.phase_final = meta.initial_rpo_phase + setup_t_nodes_phase(1) - meta.initial_epoch_time;
elseif ~isempty(setup_data) && isfield(meta, 'phase_theta_idx') && isfield(setup_data, 'y0')
    result.phase_final = setup_data.y0(meta.phase_theta_idx);
elseif isfield(meta, 'phase_seed')
    result.phase_final = meta.phase_seed;
end

result.solver = struct();
result.solver.options = [];
result.solver.stage1 = [];
result.solver.stage2 = [];
result.solver.y_stage1 = [];
result.solver.exitflag_stage1 = [];
result.solver.exitflag_stage2 = [];
result.solver.gradient_check = [];
result.solver.J_opt = [];

if ~isempty(solve_out)
    result.solver.options = solve_out.options;
    if isfield(solve_out, 'J_opt')
        result.solver.J_opt = solve_out.J_opt;
    end
    if isfield(solve_out, 'output_stage1')
        result.solver.stage1 = solve_out.output_stage1;
        result.solver.exitflag_stage1 = solve_out.exitflag_stage1;
    end
    if isfield(solve_out, 'y_stage1')
        result.solver.y_stage1 = solve_out.y_stage1;
    end
    if isfield(solve_out, 'output_stage2')
        result.solver.stage2 = solve_out.output_stage2;
        result.solver.exitflag_stage2 = solve_out.exitflag_stage2;
    end
end

if ~isempty(grad_results)
    result.solver.gradient_check = grad_results;
end

result.trajectory = [];
result.dv = [];
if ~isempty(post)
    result.trajectory = struct('t_nodes', post.t_nodes, 'xx_nodes', post.xx_nodes, 'U_nodes', post.U_nodes, ...
        'xx_traj_rot', post.xx_traj_rot, 't_traj', post.t_traj, 'STM_vec', post.STM_vec, ...
        'xx_traj_in', post.xx_traj_in, 'xx_nodes_in', post.xx_nodes_in);

    result.dv = struct('DV_mag_opt', post.DV_mag_opt, 'DVvec_opt_rot', post.DVvec_opt_rot, ...
        'DVvec_opt_in', post.DVvec_opt_in);
end

% Store the pre-optimization trajectory reconstructions used by the new plots.
result.preopt = struct();
result.preopt.setup = [];
result.preopt.stage1 = [];
if nargin >= 14 && ~isempty(preopt_results)
    if isfield(preopt_results, 'setup')
        result.preopt.setup = preopt_results.setup;
    end
    if isfield(preopt_results, 'stage1')
        result.preopt.stage1 = preopt_results.stage1;
    end
end

result.plots = plot_results;

end

function [Rp, Vp] = get_perilune_state(v_in, v_out)
%GET_PERILUNE_STATE Recover a lunar periapsis state from the two asymptotic velocities.

muM = 4902.83;

% Turning angle

dotProd = dot(v_in, v_out);
norm2 = norm(v_in)*norm(v_out);
delta = acos(dotProd/norm2);

% Pericenter radius magnitude
v_mid = (norm(v_in) + norm(v_out)) / 2;
RpMag = muM/v_mid^2 * (1/sin(delta/2) - 1);

% Perilune velocity magnitude (vis-viva)
VpMag = sqrt(v_mid^2 + 2*muM/RpMag);

% Position and velocity vectors
Rp = RpMag/(2*v_mid*sin(delta/2)) * (v_in - v_out);
Vp = VpMag/(2*v_mid*cos(delta/2)) * (v_in + v_out);
end

function [value,isterminal,direction] = soi_event_CR3BP(~, x, rSOI, mu)
%SOI_EVENT_CR3BP Stop when |r| hits the SOI radius (either direction)

posMoon = [1-mu 0 0]';
value      = norm(x(1:3) - posMoon) - rSOI;
isterminal = 1;
direction  = 0;
end

function [U, is_case2b] = dv_to_regularized_U(DV)
%DV_TO_REGULARIZED_U Convert Delta v into regularized variables U = (u,w,s)
% according to Oshima's 3D formulation.

DV = DV(:);
dvx = DV(1);
dvy = DV(2);
dvz = DV(3);

dv_mag = sqrt(dvx^2 + dvy^2 + dvz^2);
dv_yz  = dvy^2 + dvz^2;

% Very small threshold to detect zero DV

tol = 1e-14;

is_case2b = false;

if dv_mag < tol
    % Case (iii): zero DV
    fprintf('CASE~3\n');
    u = eps;
    w = eps;
    s = eps;

elseif dvx > 0
    % Case (i): dvx > 0
    fprintf('CASE~1\n');
    u = sqrt(0.5 * (dv_mag + dvx));
    w = dvy / (2*u);
    s = dvz / (2*u);

elseif (dvx <= 0) && (dv_yz > tol)
    % Case (ii-a): dvx <= 0 and (dvy^2 + dvz^2) > 0
    fprintf('CASE~2a\n');
    u = sqrt(0.5 * (dv_mag + dvx));
    w = dvy / (2*u);
    s = dvz / (2*u);

elseif (abs(dvx) > tol) && (dv_yz < tol) && (dvx < 0)
    % Case (ii-b): dvx < 0 and dvy = dvz = 0
    fprintf('CASE~2b\n');
    fprintf('CAUTION!!!\n');
    u = sqrt(-0.5 * dvx);
    w = -sqrt(-0.5 * dvx);
    s = 0.0;
    is_case2b = true;

else
    % Fallback
    fprintf('Fallback: CAUTION\n');
    u = eps;
    w = eps;
    s = eps;
end

U = [u; w; s];

end

function local_progress_log(log_config, component, message, varargin)
%LOCAL_PROGRESS_LOG Print one timestamped progress message and flush the MATLAB console.

[enabled, t_start] = local_progress_state(log_config);
if ~enabled
    return
end
time_text = char(datetime('now', 'Format', 'HH:mm:ss'));
fprintf('[%s | +%.1fs] %s | %s\n', time_text, toc(t_start), component, sprintf(message, varargin{:}));
drawnow('limitrate');
end

function [enabled, t_start] = local_progress_state(log_config)
%LOCAL_PROGRESS_STATE Resolve boolean and elapsed-time origin for console progress logs.

enabled = false;
t_start = tic;
if isstruct(log_config)
    if isfield(log_config, 'enabled') && ~isempty(log_config.enabled)
        enabled = logical(log_config.enabled);
    end
    if isfield(log_config, 't_start') && ~isempty(log_config.t_start)
        t_start = log_config.t_start;
    end
else
    enabled = ~isempty(log_config) && logical(log_config);
end
if ~isstruct(log_config) && enabled && isappdata(0, 'RPO_progress_log_t_start')
    t_start = getappdata(0, 'RPO_progress_log_t_start');
end
end
