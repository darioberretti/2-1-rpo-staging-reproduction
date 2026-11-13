function result = compute_validation(cfg,inputs)
%COMPUTE_VALIDATION Reoptimize and independently reconstruct the frozen zero-phase problem.
%
% PROTOTYPE:
%   result = sel2rpo.compute_validation(cfg, inputs)
%
% DESCRIPTION:
%   Use inputs.problems(1) without rebuilding its shooting nodes
%   from encounter databases. Reconstruct the unoptimized and solved
%   arcs separately so their actual continuity defects remain visible.
%
% INPUT:
%   cfg  Configuration with meta and sel2_to_rpo fields from
%        config/reproduction_config.m; physical units follow field suffixes.
%   inputs  Struct loaded from data/inputs/run_04_reference_inputs.mat;
%             problems(1) contains y0, nlp and constants (see
%             sel2rpo.reconstruct_transfer for the decision/state layout).
%
% OUTPUT:
%   result  Result with data.initial_guess and solution trajectories
%             from sel2rpo.reconstruct_transfer, problem, y_sol and
%             solver_diagnostics. Metrics report delta_v_m_s [m/s],
%             waiting/flight/elapsed_days [days] and the maximum absolute
%             propagated shooting defect in normalized state components.

% Reconstruct the frozen initial guess, then solve the same shooting problem.
problem = inputs.problems(1);
fprintf('Initial guess: frozen unoptimized shooting nodes, not an encounter-database rebuild.\n');
initial = sel2rpo.reconstruct_transfer(problem.y0,problem);
s = cfg.sel2_to_rpo.optimization;
case_data = sel2rpo.solve_case(problem,problem.y0,0,s,s.constraint_tolerance,1,'zero',0,1);
solution = sel2rpo.reconstruct_transfer(case_data.y_sol,problem);

% Retain both trajectories and independent propagated-defect diagnostics.
data = struct('initial_guess',initial,'solution',solution,'problem',problem, ...
    'initial_guess_kind','reconstructed_unoptimized_nodes','y_sol',case_data.y_sol, ...
    'solver_diagnostics',case_data);
metrics = struct('delta_v_m_s',solution.total_delta_v_m_s,'flight_days',solution.flight_days, ...
    'waiting_days',solution.waiting_days,'elapsed_days',solution.elapsed_days, ...
    'maximum_propagated_defect',max(abs(solution.defects),[],'all'));
result = rpo.finalize_result(struct('data',data,'metrics',metrics),cfg, ...
    'run_04_sel2_to_rpo: representative transfer','data/inputs/run_04_reference_inputs.mat: unoptimized zero problem','generated_result');
end
