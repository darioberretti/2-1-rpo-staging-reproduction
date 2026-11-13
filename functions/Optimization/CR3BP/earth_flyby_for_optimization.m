function [E, vinf_in_R, vinf_out_R] = earth_flyby_for_optimization(rp_R, vp_R, mu)
%EARTH_FLYBY_FOR_OPTIMIZATION Recover incoming/outgoing hyperbolic asymptotic velocities.
%
% PROTOTYPE:
%   [E, vinf_in_R, vinf_out_R] = earth_flyby_for_optimization(rp_R, vp_R, mu)
%
% DESCRIPTION:
%   Given an inertial periapsis state about the flyby body, use specific
%   energy and the half-deflection angle to recover the two asymptotes.
%   The caller must supply a hyperbolic periapsis (positive energy); this
%   function does not validate that physical condition.
%   The half-deflection construction follows Strange et al. (2007).
%
% INPUT:
%   rp_R  [3 elements] Body-centered inertial periapsis position [km].
%   vp_R  [3 elements] Velocity in the same inertial axes [km/s].
%   mu    Flyby-body gravitational parameter [km^3/s^2].
%
% OUTPUT:
%   E           Specific orbital energy [km^2/s^2].
%   vinf_in_R   Incoming asymptotic velocity in the supplied axes [km/s].
%   vinf_out_R  Outgoing asymptotic velocity in the supplied axes [km/s].
%                 Vector orientation follows the input vectors.

rp = norm(rp_R);     % periapsis radius
vp = norm(vp_R);     % periapsis speed

% Keplerian orbital energy at periapsis
E = 0.5 * vp^2 - mu / rp;

% Hyperbolic trajectories must satisfy E > 0
% (we don't enforce it here; caller should ensure consistency)

% Magnitude of asymptotic velocity
vinf = sqrt(2 * E);

% Half of the bending angle (delta/2)
% Strange et al. (2007) formula
delta_half = asin( mu / (mu + rp * vinf^2) );

% Local unit vectors at periapsis
e_r = rp_R / rp;      % radial unit vector
e_t = vp_R / vp;      % tangential unit vector

% Incoming v-infinity vector
vinf_in_R  = vinf * ( e_r * sin(delta_half) + e_t * cos(delta_half) );

% Outgoing v-infinity vector
vinf_out_R = vinf * (-e_r * sin(delta_half) + e_t * cos(delta_half) );

end
