function [departure_velocity_km_s, arrival_velocity_km_s, extremal_distances_km] = ...
        solve_lambert(departure_position_km, arrival_position_km, time_of_flight_days, ...
        mu_km3_s2, transfer_type)
%SOLVE_LAMBERT Solve a frozen zero-revolution Lambert problem.
%
% PROTOTYPE:
%   [v1, v2, extremal_distances_km] = rpo2ds.solve_lambert(r1, r2, ...
%       time_of_flight_days, mu_km3_s2, transfer_type)
%
% DESCRIPTION:
%   Load the bundled FEX-Lambert dependency and preserve the historical
%   direct/retrograde orientation through the sign of time of flight.
%   Require the selected zero-revolution branch to report success.
%
% INPUT:
%   departure_position_km  [3 elements] Initial position [km].
%   arrival_position_km    [3 elements] Final position [km], in the same frame.
%   time_of_flight_days    Scalar flight duration [days]; magnitude is used.
%   mu_km3_s2              Central-body gravitational parameter [km^3/s^2].
%   transfer_type          Optional direct (default) or retrograde.
%
% OUTPUT:
%   departure_velocity_km_s  [3 x 1] Departure velocity [km/s].
%   arrival_velocity_km_s    [3 x 1] Arrival velocity [km/s].
%   extremal_distances_km     Column of extremal radii from FEX-Lambert [km].

if nargin < 5 || isempty(transfer_type)
    transfer_type = "direct";
end
transfer_type = validatestring(char(transfer_type), {'direct', 'retrograde'});

persistent dependency_ready
if isempty(dependency_ready)
    repository_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    dependency_dir = fullfile(repository_root, 'third_party', 'FEX-Lambert');
    if ~isfile(fullfile(dependency_dir, 'lambert.m'))
        error('rpo2ds:MissingLambert', ...
            'FEX-Lambert v1.4 is missing from third_party/FEX-Lambert.');
    end
    addpath(dependency_dir, '-begin');
    dependency_ready = true;
end

departure_position_km = departure_position_km(:).';
arrival_position_km = arrival_position_km(:).';

% FEX-Lambert uses the sign of the time of flight to select the transfer
% orientation.  lambertMR's historical direct branch selects the same
% orientation from the z component of cross(r1,r2); preserve that choice.
cross_product = cross(departure_position_km, arrival_position_km);
orientation = sign(cross_product(3));
if orientation == 0
    orientation = 1;
end
if strcmp(transfer_type, 'retrograde')
    orientation = -orientation;
end
signed_time_of_flight_days = orientation * abs(time_of_flight_days);

[departure_velocity_km_s, arrival_velocity_km_s, extremal_distances_km, exitflag] = ...
    lambert(departure_position_km, arrival_position_km, ...
    signed_time_of_flight_days, 0, mu_km3_s2);
if exitflag ~= 1
    error('rpo2ds:LambertFailure', ...
        'FEX-Lambert failed for the frozen zero-revolution, short-way branch (exit flag %d).', ...
        exitflag);
end

departure_velocity_km_s = departure_velocity_km_s(:);
arrival_velocity_km_s = arrival_velocity_km_s(:);
extremal_distances_km = extremal_distances_km(:);
end
