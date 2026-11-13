function [v_eci, dv_dtheta] = earth_sun_rth_to_eci(v_rth, theta_sun_em, t_epoch)
%EARTH_SUN_RTH_TO_ECI Convert an Earth-Sun RTH velocity into Earth ECI.
%
% PROTOTYPE:
%   [v_eci, dv_dtheta] = earth_sun_rth_to_eci(v_rth, theta_sun_em, t_epoch)
%
% DESCRIPTION:
%   The R axis points from the Sun to Earth, T follows Earth's heliocentric
%   velocity, and H is positive out of plane. The inertial Sun angle is the
%   Earth-Moon rotating Sun phase plus the normalized epoch.
%
% INPUT:
%   v_rth         [3 elements] Velocity in the Earth-Sun RTH axes [km/s].
%   theta_sun_em  Sun angle in the Earth-Moon rotating frame [rad].
%   t_epoch       Normalized epoch [TU], numerically the frame angle [rad].
%
% OUTPUT:
%   v_eci         [3 x 1] Velocity in Earth-centered inertial axes [km/s].
%   dv_dtheta     [3 x 1] Derivative with respect to the inertial Sun angle
%                   [km/s per rad]; also the partial with respect to either
%                   supplied angle-like variable while the other is fixed.

v_rth = v_rth(:);

% Build the inertial Sun angle from the Earth-Moon rotating Sun phase.
theta_sun_eci = theta_sun_em + t_epoch;

% Build the planar Earth-Sun RTH basis in Earth-centered inertial axes.
uR = [-cos(theta_sun_eci); -sin(theta_sun_eci); 0];
uT = [ sin(theta_sun_eci); -cos(theta_sun_eci); 0];
uH = [0; 0; 1];

% Rotate the requested RTH velocity into ECI.
v_eci = v_rth(1) * uR + v_rth(2) * uT + v_rth(3) * uH;

% Build the derivative with respect to theta_sun_eci.
duR = [ sin(theta_sun_eci); -cos(theta_sun_eci); 0];
duT = [ cos(theta_sun_eci);  sin(theta_sun_eci); 0];
dv_dtheta = v_rth(1) * duR + v_rth(2) * duT;
end
