function solution = solve_frozen_flexibility_case(frozen_setup, settings)
%SOLVE_FROZEN_FLEXIBILITY_CASE Optimize one archived flexibility NLP setup.
%
% PROTOTYPE:
%   solution = rpo2ds.solve_frozen_flexibility_case(frozen_setup, settings)
%
% DESCRIPTION:
%   Build the CasADi problem from its frozen decision-vector seed and bounds,
%   run the configured feasibility/objective stages, then decode physical
%   delta-v and flight time. Compute the bound residual of the final nonlinear
%   constraints separately from the solver's reported status.
%
% INPUT:
%   frozen_setup  Struct with y0/lb/ub, problem_data, velocity_unit_km_s,
%                 and time_unit_days. Each node stores [X(6); U(3)];
%                 problem_data.dt_idx identifies segment durations [TU].
%   settings      Frozen solver iterations, tolerances, and run_stage1 flag.
%
% OUTPUT:
%   solution      Struct with converged (positive exitflag), has_solution
%                 (finite decoded payload), exitflag/objective, maximum
%                 nonlinear constraint violation, delta-v [m/s], and TOF [days].

% Solve the archived NLP with the configured iteration limits and tolerances.
opts = struct('run_stage1', settings.run_stage1, 'solver_display', 'off', ...
    'max_iter_stage1', settings.maximum_stage1_iterations, ...
    'max_iter_stage2', settings.maximum_stage2_iterations, ...
    'constraint_tolerance', settings.constraint_tolerance, ...
    'optimality_tolerance', settings.optimality_tolerance);
problem = build_casadi_rpo2ds_initial_state_problem( ...
    frozen_setup.y0, frozen_setup.lb, frozen_setup.ub, ...
    frozen_setup.problem_data, opts);
raw = solve_casadi_ipopt_two_stage(problem, frozen_setup.y0, opts);

% Decode regularized burns and segment durations into m/s and days.
y_opt = raw.y_opt(:);
node_blocks = reshape(y_opt(1:9 * frozen_setup.problem_data.n), 9, []);
delta_v_m_s = 1000 * frozen_setup.velocity_unit_km_s * ...
    sum(sum(node_blocks(7:9, :) .^ 2, 1));
time_of_flight_days = frozen_setup.time_unit_days * ...
    sum(y_opt(frozen_setup.problem_data.dt_idx));

% Evaluate final nonlinear bound violations independently of the solver status.
constraints = full(problem.constraint_function(y_opt));
lower_violation = max(problem.lbg - constraints, 0);
upper_violation = max(constraints - problem.ubg, 0);
solution = struct( ...
    'converged', isfinite(raw.exitflag) && raw.exitflag > 0, ...
    'has_solution', all(isfinite(y_opt)) && isfinite(delta_v_m_s) && ...
    isfinite(time_of_flight_days), ...
    'exitflag', raw.exitflag, ...
    'objective', raw.J_opt, ...
    'constraint_residual_inf', max([lower_violation; upper_violation]), ...
    'total_delta_v_m_s', delta_v_m_s, ...
    'time_of_flight_days', time_of_flight_days);
end
