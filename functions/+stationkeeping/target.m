function solution = target(estimate,t_start,requested_dt,target_state,cfg,sun_phase)
%TARGET Correct an impulse and coast time with a damped four-component Newton solve.
%
% PROTOTYPE:
%   solution = stationkeeping.target(estimate, t_start, requested_dt, ...
%       target_state, cfg, sun_phase)
%
% DESCRIPTION:
%   Target the frozen [x,y,z,vy] components using three initial
%   velocity increments and free coast duration. The Jacobian combines
%   STM columns mapping initial velocity to terminal state with the
%   terminal flow for time sensitivity. Apply a damped pseudoinverse
%   step and clamp duration to its requested-coast bounds. Success
%   requires every physical component tolerance, with no small-step stop.
%
% INPUT:
%   estimate      [6 x 1] Estimated EM rotating state at t_start,
%                   ordered [x;y;z;vx;vy;vz], position [DU], velocity [VU].
%   t_start       Absolute impulse epoch [TU].
%   requested_dt  Positive requested coast duration [TU].
%   target_state  [6 x 1] Fixed reference crossing state in the same units.
%   cfg           Stationkeeping settings; targeting.state_indices is
%                   [1;2;3;5]. Position tolerance [km] and vy tolerance
%                   [cm/s] are converted with cfg.units. Bounds are
%                   targeting.coast_bounds_factor*requested_dt.
%   sun_phase     Sun angle at absolute time zero [rad].
%
% OUTPUT:
%   solution      Struct with success, stop_reason, iterations,
%                   dv_plan [3 x 1, VU], coast_tu [TU], and
%                   residual_normalized [1 x 4] in [x,y,z,vy] order.
%                   Residuals are from the last evaluated iterate,
%                   divided by component tolerances; |entry|<=1
%                   satisfies that component. On failure, returned
%                   commands must not be treated as successful targets.

% Normalize physical tolerances and initialize the impulse/time correction.
indices = cfg.targeting.state_indices;
tolerances = [repmat(cfg.targeting.position_tolerance_km/cfg.units.DU_km,3,1); ...
    (cfg.targeting.vy_tolerance_cm_s*1e-5)/cfg.units.VU_km_s];
u = [0;0;0;requested_dt];
bounds = cfg.targeting.coast_bounds_factor*requested_dt;
solution = struct('success',false,'stop_reason','maximum_iterations', ...
    'iterations',0,'residual_normalized',nan(1,4),'dv_plan',zeros(3,1),'coast_tu',requested_dt);

% Propagate the commanded estimated state with an identity-initialized STM.
for iteration = 1:cfg.targeting.maximum_iterations
    corrected_state = estimate;
    corrected_state(4:6) = corrected_state(4:6)+u(1:3);
    final_time = t_start+u(4);
    out = propagate(struct('x0',corrected_state,'t0',t_start,'tf',final_time, ...
        'mu',cfg.dynamics.mu,'dimensions','3D','framework','BCR4BP', ...
        'stm_computation','true','thetaSun0',sun_phase, ...
        'absTol',cfg.propagation.absolute_tolerance,'relTol',cfg.propagation.relative_tolerance));

    % Test each targeted component and retain its tolerance-scaled residual.
    residual = out.xf-target_state;
    residual = residual(indices);
    solution.iterations = iteration;
    solution.residual_normalized = (residual./tolerances).';
    if all(abs(residual) <= tolerances)
        solution.success = true;
        solution.stop_reason = 'component_tolerances';
        break;
    end

    % Combine initial-velocity and final-time sensitivities for the Newton step.
    flow = rhs_BCR4BP_3D(final_time,out.xf,cfg.dynamics.mu,sun_phase);
    jacobian = [out.STMf(indices,4:6),flow(indices)];
    if any(~isfinite(jacobian(:))) || any(~isfinite(residual))
        solution.stop_reason = 'nonfinite_targeting';
        break;
    end

    % Damp the correction and enforce bounds relative to the requested duration.
    step = -pinv(jacobian)*residual;
    step = cfg.targeting.damping*step;
    u = u+step;
    u(4) = min(max(u(4),bounds(1)),bounds(2));
end

% Return the current command alongside its explicit convergence status.
solution.dv_plan = u(1:3);
solution.coast_tu = u(4);
end
