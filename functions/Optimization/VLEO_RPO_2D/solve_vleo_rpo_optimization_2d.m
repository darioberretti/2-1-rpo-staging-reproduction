function result = solve_vleo_rpo_optimization_2d(cfg, y0_override, seed_package, x0_rpo_planar_3bp)
%SOLVE_VLEO_RPO_OPTIMIZATION_2D Solve one planar LEO-to-RPO transfer with fmincon.
%
% PROTOTYPE:
%   result = solve_vleo_rpo_optimization_2d(cfg, y0_override, seed_package, ...
%       x0_rpo_planar_3bp)
%
% DESCRIPTION:
%   Build the mapped warm start and optionally replace it with a continuation
%   vector. Use the sensitivity-based objective/constraint pair, optionally
%   run zero-objective feasibility, then minimize regularized maneuver cost.
%   Postprocess both the initial and final vectors for geometry and residuals.
%
% INPUT:
%   cfg                   Validation solver/setup settings from cfg.earth_to_rpo.
%   y0_override           Optional continuation decision-vector column; [] uses
%                         the newly built seed. Must match the mapped layout.
%   seed_package          Optional explicit patched seed; see
%                         build_vleo_rpo_ms_initial_guess_2d.
%   x0_rpo_planar_3bp     Canonical planar target [x, y, vx, vy] in DU/VU,
%                         supplied with seed_package by the public runner.
%
% OUTPUT:
%   result                Struct with cfg/setup, initial_guess, stage1, y_opt,
%                         J_opt, exitflag_opt/output_opt, postprocessed solution,
%                         and stage runtimes [s]. Physical solution metrics use
%                         m/s, days, and km; solver flags are retained without
%                         applying a reference-agreement acceptance rule.

if nargin >= 3
    setup = build_vleo_rpo_ms_initial_guess_2d(cfg, seed_package, x0_rpo_planar_3bp);
else
    setup = build_vleo_rpo_ms_initial_guess_2d(cfg);
end

% Replace the warm start when a continuation vector is provided.
if nargin >= 2 && ~isempty(y0_override)
    setup.y0 = y0_override(:);

    % Expand the warm-start altitude bounds around continuation guesses when needed.
    if strcmpi(setup.nlp.departure_constraint_mode, 'fixed_altitude_seed_warmstart') && ~isempty(setup.nlp.map.departure_altitude_km)
        altitude_idx = setup.nlp.map.departure_altitude_km;
        setup.lb(altitude_idx) = min(setup.lb(altitude_idx), setup.y0(altitude_idx));
        setup.ub(altitude_idx) = max(setup.ub(altitude_idx), setup.y0(altitude_idx));
    end
end

% Build the fast objective and constraint handles.
objfun_handle = @(y) objfun_vleo_rpo_regularized_2d_fast(y, setup.nlp);
constr_handle = @(y) constr_vleo_rpo_regularized_2d_fast(y, setup.nlp);

% Postprocess the active initial guess for diagnostics.
initial_guess = postprocess_vleo_rpo_optimization_2d(setup.y0, setup, 'initial_guess');

% Build the stage options once for the two fmincon calls.
options_stage1 = local_build_fmincon_options(cfg.solver_stage1, cfg.max_iter_stage1, cfg.max_fun_evals_stage1, cfg.optimality_tolerance_stage1, cfg);
options_stage2 = local_build_fmincon_options(cfg.solver_stage2, cfg.max_iter_stage2, cfg.max_fun_evals_stage2, cfg.optimality_tolerance_stage2, cfg);

% Initialize the stage-1 container before the feasibility solve.
stage1 = struct;
stage1.y = setup.y0;
stage1.fval = 0;
stage1.exitflag = NaN;
stage1.output = struct;
stage1.runtime_s = 0;

% Run the zero-objective feasibility stage when requested.
if cfg.use_zero_objective_stage
    stage1_tic = tic;
    [stage1.y, stage1.fval, stage1.exitflag, stage1.output] = fmincon(@(y) obj_zero(y), setup.y0, [], [], [], [], setup.lb, setup.ub, constr_handle, options_stage1);
    stage1.runtime_s = toc(stage1_tic);
end

% Run the final delta-v optimization stage from the stage-1 output.
stage2_tic = tic;
[y_opt, J_opt, exitflag_opt, output_opt] = fmincon(objfun_handle, stage1.y, [], [], [], [], setup.lb, setup.ub, constr_handle, options_stage2);
stage2_runtime_s = toc(stage2_tic);

% Reconstruct and summarize the optimized transfer.
solution = postprocess_vleo_rpo_optimization_2d(y_opt, setup, 'optimized_solution');

% Package the complete one-altitude solve result.
result = struct;
result.cfg = cfg;
result.setup = setup;
result.initial_guess = initial_guess;
result.stage1 = stage1;
result.y_opt = y_opt;
result.J_opt = J_opt;
result.exitflag_opt = exitflag_opt;
result.output_opt = output_opt;
result.stage2_runtime_s = stage2_runtime_s;
result.solution = solution;
result.runtime_s = stage1.runtime_s + stage2_runtime_s;

end

function options = local_build_fmincon_options(algorithm_name, max_iter, max_fun_evals, optimality_tolerance, cfg)
%LOCAL_BUILD_FMINCON_OPTIONS Build the configured derivative-aware fmincon stage options.
%
% PROTOTYPE:
%   options = local_build_fmincon_options(algorithm_name, max_iter, ...
%       max_fun_evals, optimality_tolerance, cfg)
%
% DESCRIPTION:
%   Require supplied objective/constraint gradients and retain the shared
%   finite-difference, display, step, and constraint-tolerance settings.
%
% INPUT:
%   algorithm_name        fmincon algorithm for the selected stage.
%   max_iter/max_fun_evals Stage iteration/evaluation limits.
%   optimality_tolerance  Stage-specific optimality tolerance.
%   cfg                   Shared display, constraint, and step settings.
%
% OUTPUT:
%   options               optimoptions object passed to one fmincon call.

options = optimoptions('fmincon', ...
    'Algorithm', algorithm_name, ...
    'Display', cfg.display, ...
    'SpecifyObjectiveGradient', true, ...
    'SpecifyConstraintGradient', true, ...
    'FiniteDifferenceType', 'central', ...
    'MaxIterations', max_iter, ...
    'MaxFunctionEvaluations', max_fun_evals, ...
    'ConstraintTolerance', cfg.constraint_tolerance, ...
    'StepTolerance', cfg.step_tolerance, ...
    'OptimalityTolerance', optimality_tolerance);

end
