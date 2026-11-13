function costs = build_pair(radius_au, angle_deg, settings, constants)
%BUILD_PAIR Evaluate cumulative Lambert excess-speed requirements at one geometry.
%
% PROTOTYPE:
%   costs = comet.build_pair(radius_au, angle_deg, settings, constants)
%
% DESCRIPTION:
%   Transform rotating departure/intercept positions to inertial states.
%   At each flight time, inspect branches in revolution order and stop at
%   the first unavailable branch. Exclude negative semimajor axes and retain
%   strict cost improvements over the available flight-time interval.
%
% INPUT:
%   radius_au  Intercept radius [AU].
%   angle_deg  Intercept angle in rotating heliocentric x-y axes [deg].
%   settings   tof_grid_days [T values, days], max_revolutions, and
%              departure_rot_au [3 x 1, AU].
%   constants  AU_km, year_days, year_s, n_earth_rads and mu_sun_km3s2.
%
% OUTPUT:
%   costs      [T x 1 double] Minimum departure excess speed [km/s] over
%              flight times up to each grid value. Unavailable costs are NaN.

times = settings.tof_grid_days(:);
costs = NaN(numel(times), 1);
cumulative_best = Inf;

% Set the corotating departure reference and fixed rotating intercept.
[r1, reference_velocity] = comet_rotating_frame_state(settings.departure_rot_au, 0, constants);
arrival_rot = [radius_au * cosd(angle_deg); radius_au * sind(angle_deg); 0];

for i_time = 1:numel(times)
    [r2, ~] = comet_rotating_frame_state(arrival_rot, ...
        times(i_time) / constants.year_days, constants);
    stop = false;
    best = Inf;

    % Preserve branch availability and elliptic-conic eligibility.
    for rev = 0:settings.max_revolutions
        if stop, break; end
        for branch = 0:(rev > 0)
            [velocity, ~, a, ~, flag] = comet.lambert_branch( ...
                r1, r2, times(i_time), constants.mu_sun_km3s2, rev, branch);
            if flag ~= 1, stop = true; break; end
            if a < 0, continue; end
            excess = velocity - reference_velocity;
            cost = norm(excess);
            if cost < best, best = cost; end
        end
    end

    % Equal costs retain the earlier flight time's minimum.
    if best < cumulative_best, cumulative_best = best; end
    if cumulative_best < Inf, costs(i_time) = cumulative_best; end
end
end
