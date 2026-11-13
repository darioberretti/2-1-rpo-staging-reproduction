function [hit, time, state, trajectory_time, trajectory] = propagate_boundary(magnitude, angle, altitude_km, cfg, rpo_state)
%PROPAGATE_BOUNDARY Propagate a planar Earth-access candidate backward to its boundary.
%
% PROTOTYPE:
%   [hit, time, state, trajectory_time, trajectory] = ...
%       earth2rpo.propagate_boundary(magnitude, angle, altitude_km, cfg, rpo_state)
%
% DESCRIPTION:
%   Add the selected impulse to the RPO state and integrate the BCR4BP with
%   initial Sun phase zero. Stop at the requested Earth altitude or either
%   body-clearance limit. The target-altitude event is active only for
%   positive Earth-centered inertial angular momentum.
%
% INPUT:
%   magnitude    Scalar impulse magnitude in normalized VU.
%   angle        Impulse direction in the rotating plane [rad].
%   altitude_km  Target Earth-boundary altitude [km].
%   cfg          cfg.earth_to_rpo; shared supplies historical DU/VU/TU scales
%                  and body sizes, access supplies horizon/tolerances/clearances.
%   rpo_state    [4 elements] Rotating barycentric [x, y, vx, vy] in DU/VU.
%
% OUTPUT:
%   hit              True when the terminal event is the requested boundary.
%   time             Signed backward event time [TU], or NaN on a miss.
%   state            [1 x 4] Terminal state in DU/VU, or NaNs on a miss.
%   trajectory_time  Full ode113 time column [TU], decreasing from zero.
%   trajectory       Matching [Nt x 4] rotating states in DU/VU.

s = cfg.shared;
a = cfg.access;

% Apply the departure impulse in rotating axes before backward integration.
initial_state = rpo_state(:);
initial_state(3:4) = initial_state(3:4) + magnitude * [cos(angle); sin(angle)];

% Integrate until the target boundary, a body-clearance event, or the time horizon.
options = odeset('RelTol', a.relative_tolerance, 'AbsTol', a.absolute_tolerance, ...
    'Events', @(t,y) boundary_events(t, y, s, a, altitude_km));
[trajectory_time, trajectory, event_time, event_state, event_index] = ode113( ...
    @(t,y) rhs_BCR4BP_2D(t,y,s.mu,0), [0, -a.horizon_days/s.TU_days], initial_state, options);

% Accept only the target event; keep full trajectories available even on a miss.
hit = ~isempty(event_index) && event_index(end) == 3;
time = NaN;
state = NaN(1,4);
if hit
    time = event_time(end);
    state = event_state(end,:);
end
end

function [value, terminal, direction] = boundary_events(~, y, s, a, altitude_km)
%BOUNDARY_EVENTS Identify Earth/Moon clearance violations and the requested backward boundary.
%
% PROTOTYPE:
%   [value, terminal, direction] = boundary_events(t, y, s, a, altitude_km)
%
% DESCRIPTION:
%   Keep the target event nonzero unless the Earth-centered inertial angular
%   momentum is positive. All three events terminate propagation in either
%   crossing direction; the caller accepts only event index 3.
%
% INPUT:
%   t           Unused current integration time.
%   y           [4 x 1] Rotating barycentric state in DU/VU.
%   s           cfg.shared scales, mass ratio, and body radii.
%   a           cfg.access clearance settings [km].
%   altitude_km Requested Earth-boundary altitude [km].
%
% OUTPUT:
%   value       [3 x 1] Moon clearance, Earth clearance, target-altitude
%               residuals in DU; inactive target value is 1.
%   terminal    [1; 1; 1].
%   direction   [0; 0; 0].

earth_distance = norm(y(1:2) - [-s.mu;0]);
moon_distance = norm(y(1:2) - [1-s.mu;0]);
angular_momentum = (y(1)+s.mu)*(y(4)+y(1)+s.mu) - y(2)*(y(3)-y(2));
target_event = 1;
if angular_momentum > 0
    target_event = earth_distance - (s.R_e/s.DU + altitude_km/s.DU);
end
value = [moon_distance - (s.R_m/s.DU+a.moon_clearance_km/s.DU); ...
    earth_distance - (s.R_e/s.DU+a.earth_clearance_km/s.DU); target_event];
terminal = [1;1;1];
direction = [0;0;0];
end
