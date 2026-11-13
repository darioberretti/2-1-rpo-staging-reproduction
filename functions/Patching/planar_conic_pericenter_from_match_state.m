function [rp_vec, vp_vec, conic] = planar_conic_pericenter_from_match_state(r_match, v_match, muE, n_samples)
%PLANAR_CONIC_PERICENTER_FROM_MATCH_STATE Recover periapsis and the physical arc of an Earth-centered planar conic.
%
% PROTOTYPE:
%   [rp_vec, vp_vec, conic] = planar_conic_pericenter_from_match_state( ...
%       r_match, v_match, muE, n_samples)
%
% DESCRIPTION:
%   Recover energy, angular momentum, and eccentricity from the match state,
%   then construct the nondegenerate conic periapsis and its match-point time.
%   Use elliptic, hyperbolic, or parabolic-limit anomaly formulas. Periapsis
%   altitude uses the historical astroConstants(23) Earth radius, rather than
%   silently substituting the separate transfer-configuration radius.
%
% INPUT:
%   r_match   [2 elements] Earth-centered inertial match position [km].
%   v_match   [2 elements] Velocity in the same inertial axes [km/s].
%   muE       Earth gravitational parameter [km^3/s^2].
%   n_samples Optional arc sample count, default 400; values below 2 skip
%               propagation and retain empty sampled arrays.
%
% OUTPUT:
%   rp_vec, vp_vec  [1 x 2] Inertial periapsis position [km] and velocity [km/s].
%   conic           Struct with energy [km^2/s^2], h [km^2/s], eccentricity,
%                     semimajor axis/semilatus rectum/radii/altitude [km],
%                     true anomaly [rad], periapsis-to-match time [s], and
%                     optional t_arc [s], rr_arc [km], vv_arc [km/s].
%                     Sampled endpoint errors use km and km/s.

if nargin < 4
    n_samples = 400;
end

% Store Earth radius for the periapsis altitude.
R_e = astroConstants(23);

% Force row vectors for planar algebra.
r_match = reshape(r_match, 1, []);
v_match = reshape(v_match, 1, []);

% Compute the basic scalar invariants.
r_norm = norm(r_match);
v_norm = norm(v_match);
rv = dot(r_match, v_match);
h = r_match(1) * v_match(2) - r_match(2) * v_match(1);
energy = 0.5 * v_norm^2 - muE / r_norm;

% Compute the eccentricity vector and its magnitude.
e_vec = ((v_norm^2 - muE / r_norm) * r_match - rv * v_match) / muE;
ecc = norm(e_vec);

% Recover the main conic parameters.
p = h^2 / muE;
a = -muE / (2 * energy);
rp_radius = p / (1 + ecc);
rp_altitude = rp_radius - R_e;

% Recover the periapsis direction from the eccentricity vector.
rp_hat = e_vec / ecc;
t_hat = sign(h) * [-rp_hat(2), rp_hat(1)];

% Recover the periapsis state.
vp_mag = sqrt(2 * (energy + muE / rp_radius));
rp_vec = rp_radius * rp_hat;
vp_vec = vp_mag * t_hat;

% Recover the match-point true anomaly.
cos_nu = dot(e_vec, r_match) / (ecc * r_norm);
sin_nu = rv * h / (muE * ecc * r_norm);
true_anomaly_match = atan2(sin_nu, cos_nu);

% Recover the time from periapsis to the boundary point.
tof_peri_to_match = local_time_from_true_anomaly(true_anomaly_match, ecc, a, p, muE);

% Sample the physical branch only when the caller requests it.
tt_arc = [];
rr_arc = [];
vv_arc = [];
match_position_error = NaN;
match_velocity_error = NaN;

if n_samples >= 2
    t_arc = linspace(0, tof_peri_to_match, n_samples).';
    x0 = [rp_vec(:); vp_vec(:)];
    [tt_arc, xx_arc] = ode113(@(t, x) local_rhs_2bp_planar(t, x, muE), t_arc, x0);
    rr_arc = xx_arc(:, 1:2);
    vv_arc = xx_arc(:, 3:4);
    match_position_error = norm(rr_arc(end, :) - r_match);
    match_velocity_error = norm(vv_arc(end, :) - v_match);
end

% Store the output structure.
conic = struct;
conic.energy = energy;
conic.h = h;
conic.eccentricity = ecc;
conic.eccentricity_vector = e_vec;
conic.semimajor_axis = a;
conic.semilatus_rectum = p;
conic.rp_radius = rp_radius;
conic.rp_altitude = rp_altitude;
conic.true_anomaly_match = true_anomaly_match;
conic.tof_peri_to_match = tof_peri_to_match;
conic.t_arc = tt_arc;
conic.rr_arc = rr_arc;
conic.vv_arc = vv_arc;
conic.match_position_error = match_position_error;
conic.match_velocity_error = match_velocity_error;
conic.radial_velocity_match = rv / r_norm;
conic.is_hyperbolic = ecc > 1;
conic.is_elliptic = ecc < 1;
conic.is_parabolic = abs(ecc - 1) <= 1e-10;

end

function dx = local_rhs_2bp_planar(~, x, muE)
%LOCAL_RHS_2BP_PLANAR Evaluate planar two-body dynamics for [r(2); v(2)] in km and seconds.

r = x(1:2);
v = x(3:4);
r_norm = norm(r);

dx = [v; -muE * r / r_norm^3];

end

function tof = local_time_from_true_anomaly(nu, ecc, a, p, muE)
%LOCAL_TIME_FROM_TRUE_ANOMALY Convert match anomaly to periapsis-to-match time.
%
% DESCRIPTION:
%   nu is in rad, eccentricity is dimensionless, a/p are in km, and
% muE is in km^3/s^2. Return time in seconds using the active conic branch.

if ecc < 1 - 1e-10
    % Use the elliptic anomaly for bound conics.
    E = 2 * atan2(sqrt(1 - ecc) * sin(nu / 2), sqrt(1 + ecc) * cos(nu / 2));
    M = E - ecc * sin(E);
    if M < 0
        M = M + 2 * pi;
    end
    n = sqrt(muE / a^3);
    tof = M / n;
elseif ecc > 1 + 1e-10
    % Use the hyperbolic anomaly for escape conics.
    F = 2 * atanh(sqrt((ecc - 1) / (ecc + 1)) * tan(nu / 2));
    Mh = ecc * sinh(F) - F;
    n = sqrt(muE / (-a)^3);
    tof = Mh / n;
else
    % Use Barker's equation for the parabolic limit.
    D = tan(nu / 2);
    tof = 0.5 * sqrt(p^3 / muE) * (D + D^3 / 3);
end

end
