function [r_in_km, v_in_kms] = comet_rotating_frame_state(r_rot_au, epoch_years, constants)
%COMET_ROTATING_FRAME_STATE Map a fixed rotating heliocentric point to an inertial corotating state.
%
% PROTOTYPE:
%   [r_in_km, v_in_kms] = comet_rotating_frame_state( ...
%       r_rot_au, epoch_years, constants)
%
% DESCRIPTION:
%   Rotate about +z at the configured Earth angular rate from
%   aligned axes at epoch zero. Return the velocity of the rotating
%   reference point, n cross r, for excess-speed calculations; this
%   does not propagate a free Keplerian spacecraft trajectory.
%
% INPUT:
%   r_rot_au    [3 elements] Fixed position in the uniformly rotating
%                 heliocentric frame [AU].
%   epoch_years Scalar epoch [years] from aligned axes.
%   constants   Struct with n_earth_rads [rad/s], year_s [s/year]
%                 and AU_km [km/AU].
%
% OUTPUT:
%   r_in_km     [3 x 1] Heliocentric inertial position [km].
%   v_in_kms    [3 x 1] Inertial corotating reference velocity [km/s].

    phi_rad = constants.n_earth_rads * epoch_years * constants.year_s; % inertial frame rotation angle [rad]
    rotation = [cos(phi_rad), -sin(phi_rad), 0.0; ...
        sin(phi_rad),  cos(phi_rad), 0.0; ...
        0.0,           0.0,          1.0];
    r_in_au = rotation * r_rot_au(:);                                  % inertial position [AU]
    r_in_km = constants.AU_km * r_in_au;                               % inertial position [km]
    v_in_kms = constants.n_earth_rads * [-r_in_km(2); r_in_km(1); 0.0]; % inertial corotating reference velocity [km/s]
end
