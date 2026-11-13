function derivative = ode_2bp(~, state, mu)
%ODE_2BP Evaluate two-body Cartesian position/velocity derivatives.
%
% PROTOTYPE:
%   derivative = ode_2bp(t, state, mu)
%
% DESCRIPTION:
%   Use consistent length/time units; Run02 supplies km and seconds.
%
% INPUT:
%   t      Time input, unused by the autonomous equations.
%   state  [6 x 1] [r(3); v(3)], using km and km/s in Run02.
%   mu     Central-body GM [km^3/s^2] for those physical units.
%
% OUTPUT:
%   derivative  [6 x 1] [v; acceleration], with acceleration in km/s^2.

position = state(1:3);
velocity = state(4:6);
derivative = [velocity; -mu * position / norm(position)^3];
end
