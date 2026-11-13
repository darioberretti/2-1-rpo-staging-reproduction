function case_data = solve_case(problem,y_initial,theta,settings,tolerance,problem_index,family,direction,index)
%SOLVE_CASE Run feasibility SQP followed by interior-point delta-v minimization.
%
% PROTOTYPE:
%   case_data = sel2rpo.solve_case(problem, y_initial, theta, settings, ...
%       tolerance, problem_index, family, direction, index)
%
% DESCRIPTION:
%   Use exact supplied objective/constraint gradients with the
%   historical two-stage fmincon settings. Keep signed/unwrapped Sun
%   phases for families 1/2 and wrap family 3. Return solver flags and
%   residuals without rejecting a case on their values.
%
% INPUT:
%   problem  Frozen shooting problem from inputs.problems in
%            data/inputs/run_04_reference_inputs.mat. Contains nlp
%            (node maps, boundary seeds, phases, model constants,
%            constraints and ODE tolerances) and constants
%            (TU_days [days/TU_EM], VU_km_s [km/s per VU_EM]).
%   y_initial  Initial decision vector; layout in sel2rpo.reconstruct_transfer.
%   theta      Sun phase at EM time zero [rad].
%   settings   cfg.sel2_to_rpo.optimization settings.
%   tolerance  Constraint tolerance passed to both fmincon stages.
%   problem_index  Seed-family index (1, 2 or 3); 3 enables phase wrapping.
%   family     Text label for logs and result provenance.
%   direction  Branch label: -1, 0 (center) or +1.
%   index      Case index within the branch.
%
% OUTPUT:
%   case_data  Struct with y_sol, problem/family/branch labels, actual
%                solver phase [rad], delta_v_m_s [m/s], waiting/flight/
%                elapsed_days [days], maximum equality/inequality
%                residuals in normalized constraints, both exitflags
%                and solver_outputs, tolerance and runtime_seconds [s].

% Require the optimizer only for fresh transfer solves.
if isempty(which('fmincon'))
    error('sel2rpo:OptimizationToolbox','Fresh Run 04 transfers require Optimization Toolbox (fmincon). Saved-result replay does not.');
end
% The original zero/pi campaign used signed/unwrapped phases. Only the later
% 1.3-rad campaign and its extension wrapped phases before optimization.
nlp = problem.nlp; nlp.theta_sun0_rad = theta;
if problem_index==3, nlp.theta_sun0_rad=mod(theta,2*pi); end
constraints = @(y) constr_SEL2RPO_regularized_phase(y,nlp);
objective = @(y) objfun_SEL2RPO_regularized_phase(y,nlp);

% Configure the shared tolerances before the two sequential solver stages.
options = optimoptions('fmincon','SpecifyObjectiveGradient',true, ...
    'SpecifyConstraintGradient',true,'FiniteDifferenceType','central', ...
    'MaxIterations',settings.maximum_iterations,'MaxFunctionEvaluations',settings.maximum_evaluations, ...
    'OptimalityTolerance',settings.optimality_tolerance,'StepTolerance',settings.step_tolerance, ...
    'ConstraintTolerance',tolerance,'Display',settings.display);
clock_start = tic;
fprintf('%s: Sun phase %.9f rad, branch %+d, case %d; feasibility stage\n',family,nlp.theta_sun0_rad,direction,index);
options.Algorithm = 'sqp';
[y_feasible,~,flag1,output1] = fmincon(@obj_zero,y_initial,[],[],[],[],[],[],constraints,options);
fprintf('%s: Sun phase %.9f rad; DV minimization stage\n',family,nlp.theta_sun0_rad);
options.Algorithm = 'interior-point';
[y,~,flag2,output2] = fmincon(objective,y_feasible,[],[],[],[],[],[],constraints,options);

% Evaluate final constraint residuals and decode physical cost/time metrics.
[c,ceq] = constraints(y);
waiting = mod(y(nlp.map.theta_sel2)-nlp.theta_halo0_nd,nlp.sel2_period_nd)/nlp.omega_se_rad_s/86400;
flight = sum(y(nlp.map.dt_nodes))*problem.constants.TU_days;
impulses = zeros(nlp.N,3);
for node=1:nlp.N
    u = y(nlp.map.node(node).U);
    impulses(node,:) = [u(1)^2-u(2)^2-u(3)^2,2*u(1)*u(2),2*u(1)*u(3)];
end
case_data = struct('problem_index',problem_index,'theta_sun0_rad',nlp.theta_sun0_rad, ...
    'solver_theta_sun0_rad',nlp.theta_sun0_rad, ...
    'y_sol',y,'delta_v_m_s',sum(vecnorm(impulses,2,2))*problem.constants.VU_km_s*1000, ...
    'flight_days',flight,'waiting_days',waiting,'elapsed_days',waiting+flight, ...
    'equality_residual',max(abs(ceq)),'inequality_violation',max(max(c,0)), ...
    'exitflags',[flag1,flag2],'source_family',family,'branch_direction',direction, ...
    'source_case_index',index,'constraint_tolerance',tolerance, ...
    'solver_outputs',{{output1,output2}},'runtime_seconds',toc(clock_start));
fprintf('%s: phase %.9f, flags %d/%d, DV %.6f m/s, equality %.3e, inequality %.3e; %.1f s\n', ...
    family,nlp.theta_sun0_rad,flag1,flag2,case_data.delta_v_m_s, ...
    case_data.equality_residual,case_data.inequality_violation,case_data.runtime_seconds);
end
