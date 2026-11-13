function [x_rot, dx_dphase, dx_dh] = circular_parking_state_2d(phase_departure, altitude_km, constants)
%CIRCULAR_PARKING_STATE_2D Build a prograde planar circular departure state and derivatives.
%
% PROTOTYPE:
%   [x_rot, dx_dphase, dx_dh] = circular_parking_state_2d( ...
%       phase_departure, altitude_km, constants)
%
% DESCRIPTION:
%   Construct the Earth-centered inertial parking circle and transform its
%   state into rotating barycentric coordinates at epoch zero. Differentiate
%   the boundary state analytically with respect to phase and physical altitude.
%
% INPUT:
%   phase_departure  Circular-orbit position angle [rad].
%   altitude_km      Altitude above constants.R_e [km].
%   constants        Struct with R_e/DU [km], VU [km/s], muE [km^3/s^2],
%                      and Earth-Moon mass ratio mu [-].
%
% OUTPUT:
%   x_rot            [4 x 1] [x; y; vx; vy] in rotating DU/VU.
%   dx_dphase        [4 x 1] Boundary-state derivative per radian.
%   dx_dh            [4 x 1] Boundary-state derivative per km of altitude.

r_km = constants.R_e + altitude_km;
r_du = r_km / constants.DU;

v_circ_kms = sqrt(constants.muE / r_km);
v_circ_vu = v_circ_kms / constants.VU;

cphi = cos(phase_departure);
sphi = sin(phase_departure);

x_in = r_du * cphi;
y_in = r_du * sphi;

vx_in = -v_circ_vu * sphi;
vy_in = v_circ_vu * cphi;

% Convert the departure state into the rotating frame at t = 0.
x_rot = [x_in - constants.mu; y_in; vx_in + y_in; vy_in - x_in];

% Differentiate the rotating state with respect to the phase.
dx_dphase = [-r_du * sphi; r_du * cphi; cphi * (r_du - v_circ_vu); -sphi * (v_circ_vu - r_du)];

% Differentiate the circular speed with respect to altitude.
dv_dh_kms = -0.5 * v_circ_kms / r_km;
dv_dh_vu = dv_dh_kms / constants.VU;

% Differentiate the rotating state with respect to altitude.
dx_dh = [cphi / constants.DU; sphi / constants.DU; sphi * (1 / constants.DU - dv_dh_vu); cphi * (dv_dh_vu - 1 / constants.DU)];

end
