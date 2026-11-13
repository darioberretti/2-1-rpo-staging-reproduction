function result = compute_validation(cfg, map, inputs)
%COMPUTE_VALIDATION Build the patch and solve one fixed-altitude LEO-to-RPO transfer.
%
% PROTOTYPE:
%   result = earth2rpo.compute_validation(cfg, map, inputs)
%
% DESCRIPTION:
%   Require fmincon, reconstruct the historical patch, and solve the selected
%   validation altitude from a fresh seed. Record both solver stages and
%   postprocessed metrics; this function does not run an altitude continuation
%   or apply an independent reference-agreement acceptance test.
%
% INPUT:
%   cfg     cfg.earth_to_rpo, including the direct validation altitude/options.
%   map     Selected 10000 km Earth-access candidate table.
%   inputs  Tracked x0_rpo_bcr4bp and x0_rpo_3bp, each [4 x 1] in DU/VU.
%
% OUTPUT:
%   result  Struct with seed, initial_setup/initial_guess, setup, solution,
%             stage1, final exit flag/output under data; metrics contain
%             delta-v [m/s], TOF [days], sampled lunar altitude [km], and
%             solver runtime [s]. Configuration/provenance are added later.

% Require the solver before reconstructing the seed and starting optimization.
if ~license('test','Optimization_Toolbox') || isempty(which('fmincon'))
    error('earth2rpo:OptimizationToolbox', ...
        'Fresh transfer optimization requires Optimization Toolbox. Saved-result replay does not.');
end

% Build the patch and solve this altitude without a continuation override.
seed = earth2rpo.build_seed(map,cfg,inputs);
solve = solve_vleo_rpo_optimization_2d(cfg.validation,[],seed,inputs.x0_rpo_3bp);

% Retain the full setup/stage payload separately from physical summary metrics.
result.data = struct('seed',seed,'initial_setup',solve.setup,'initial_guess',solve.initial_guess, ...
    'setup',solve.setup,'solution',solve.solution,'stage1',solve.stage1, ...
    'exitflag_stage2',solve.exitflag_opt,'solver_output',solve.output_opt);
result.metrics = struct('total_dv_mps',solve.solution.dv_total_mps, ...
    'departure_dv_mps',solve.solution.departure_summary.impulse_magnitude_mps, ...
    'insertion_dv_mps',solve.solution.arrival_summary.impulse_magnitude_mps, ...
    'tof_days',solve.solution.tof_days,'minimum_moon_altitude_km',solve.solution.min_moon_altitude_km, ...
    'runtime_s',solve.runtime_s);
end
