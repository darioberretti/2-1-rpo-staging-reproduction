function result = compute_validation(cfg, database, seed, target_name, target)
%COMPUTE_VALIDATION Run one frozen Venus/Mars fixed-target validation chain.
%
% PROTOTYPE:
%   result = rpo2ds.compute_validation(cfg, database, seed, target_name, target)
%
% DESCRIPTION:
%   Construct the patched seed and solve the configured regularized transfer
%   with the phase-aware database. Stop if the solver exit status is not
%   positive; numerical agreement with a reference is assessed separately.
%
% INPUT:
%   cfg          Configuration from reproduction_config(); uses
%                rpo_to_deep_space.validation.shared and the target profile.
%   database     Phase-aware RPO-Moon table from rpo2ds.load_database.
%   seed         Struct with state_2d = [x; y; vx; vy] in rotating DU/VU.
%   target_name  venus or mars.
%   target       Corresponding target struct from rpo2ds.load_targets;
%                Earth-relative ECI/RTH v-infinity vectors are in km/s.
%
% OUTPUT:
%   result       Compact transfer result: data follows rpo2ds.pack_transfer_result,
%                with metrics, convergence, status, and producer provenance.
%                Configuration/input metadata are attached by finalize_result.

% Read the frozen target profile and require its numerical backend.
settings = cfg.rpo_to_deep_space.validation.shared;
case_settings = cfg.rpo_to_deep_space.validation.(target_name);
overall_timer = tic;
rpo2ds.require_solver(settings.solver_backend);

% Build the explicit geometry, timing, target, and solve options.
params = rpo2ds.make_transfer_params(settings, target, 'eci', 0);
params.maxTOF = case_settings.maximum_total_time_days;
params.maxTransferTOF = settings.maximum_transfer_time_days;
params.operational_mode = 'standard';
params.phase_mode = 'no_phase';
params.target_frame = 'eci';
params.vinf_target_rth_km_s = target.vinf_earth_rth_km_s(:);

opts = struct();
opts.graph = struct('flyby', false, 'initial_guess', false, 'nodes', false, ...
    'initial_guess_fmincon_rot', false, 'initial_guess_fmincon_inertial', false, ...
    'optimized_rot', false, 'optimized_inertial', false, ...
    'objzero_rot', false, 'objzero_inertial', false);
opts.solver_backend = settings.solver_backend;
opts.run_stage1 = settings.run_stage1;
opts.max_iter_stage1 = settings.maximum_stage1_iterations;
opts.max_iter_stage2 = settings.maximum_stage2_iterations;
opts.constraint_tolerance = settings.constraint_tolerance;
opts.optimality_tolerance = settings.optimality_tolerance;
opts.print_logs = true;
opts.print_patch_logs = false;
opts.progress_logs = false;
opts.setup_only = false;
opts.solver_display = 'off';

% Construct the patched seed and solve the configured two-stage NLP.
fprintf('Running %s validation with %s and the phase-aware database.\n', ...
    upper(target_name), settings.solver_backend);
raw = RPO2DS_transfer(settings.solve_mode, target.vinf_earth_eci_km_s(:), ...
    database, seed.state_2d(:), params, opts);

% Keep the reusable numerical payload and require a positive solver exit flag.
packed = rpo2ds.pack_transfer_result(raw, target_name, target, settings);
if ~packed.convergence.completed
    error('rpo2ds:ValidationFailed', ...
        '%s validation did not converge (exit flag %g). Review the solver output before rerunning.', ...
        upper(target_name), packed.convergence.exitflag);
end

% Attach section metrics and solver status before shared provenance is added.
result = struct();
result.data = packed;
result.metrics = packed.metrics;
result.convergence = packed.convergence;
result.status = struct('completed', packed.convergence.completed, ...
    'failed_sections', strings(0, 1));
result.provenance.historical_producer = "RPO2DS_transfer / frozen validation profile";
result.provenance.elapsed_seconds = toc(overall_timer);
end
