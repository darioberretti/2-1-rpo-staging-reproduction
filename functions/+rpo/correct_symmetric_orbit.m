function orbit = correct_symmetric_orbit(x0, vy0, crossing_number, cfg, mu)
%CORRECT_SYMMETRIC_ORBIT Correct vy0 so vx vanishes at a y = 0 crossing.
%
% PROTOTYPE:
%   orbit = rpo.correct_symmetric_orbit(x0, vy0, crossing_number, cfg, mu)
%
% DESCRIPTION:
%   Hold x0 fixed and use the event-time-corrected STM sensitivity to update
%   vy0. Chain segment STMs from the initial state to the selected crossing.
%   States use the normalized Earth-Moon rotating frame.
%
% INPUT:
%   x0               Scalar initial x [DU].
%   vy0              Scalar initial y-velocity guess [VU].
%   crossing_number  Positive integer selecting the terminal y = 0 crossing.
%   cfg              reference_rpo.a1 settings from reproduction_config();
%                      supplies iteration, residual, step, event, and time limits.
%   mu               Secondary-to-total mass ratio [-].
%
% OUTPUT:
%   orbit            Struct with valid, state [4 x 1] = [x; y; vx; vy],
%                      t_half [TU], x_final [DU], iterations, and residual [VU].
%                      valid remains false if a crossing fails or correction
%                      reaches its iteration limit.

state = [x0; 0; 0; vy0];
previous_vy = 0;
orbit = struct('valid', false, 'state', state, 't_half', NaN, ...
    'x_final', NaN, 'iterations', 0, 'residual', Inf);

% Correct the initial velocity using the selected crossing and its STM sensitivity.
for iteration = 1:cfg.maximum_iterations
    [crossings, times, valid] = rpo.propagate_crossings(state, ...
        crossing_number, cfg.correction_time_limit_tu, mu, true, ...
        cfg.event_tolerance);
    if ~valid
        return
    end
    final = crossings{end}.state;

    % Compose segment STMs to map the initial perturbation to the final crossing.
    stm = eye(4);
    for crossing = 1:crossing_number
        stm = crossings{crossing}.stm * stm;
    end
    derivative = rhs_CR3BP_2D(times(end), final, mu);

    % Include the crossing-time sensitivity so the Newton step stays on y = 0.
    correction_derivative = stm(3, 4) ...
        - derivative(3) / derivative(2) * stm(2, 4);
    residual = -final(3);
    step = state(4) - previous_vy;
    previous_vy = state(4);
    state(4) = state(4) + correction_derivative \ residual;

    % Accept only when both the crossing residual and previous velocity step converge.
    if abs(residual) <= cfg.residual_tolerance ...
            && abs(step) <= cfg.step_tolerance
        orbit.valid = true;
        orbit.state = state;
        orbit.t_half = times(end);
        orbit.x_final = final(1);
        orbit.iterations = iteration;
        orbit.residual = abs(residual);
        return
    end
end
end
