function [rows, retried] = compute_task(axis, departure_state, settings)
%COMPUTE_TASK Propagate one impulsive departure and retain every lunar entry.
%
% PROTOTYPE:
%   [rows, retried] = rpomoon.compute_task(axis, departure_state, settings)
%
% DESCRIPTION:
%   Preserve the historical forward/backward CR3BP encounter calculation.
%   Retry integration once at relaxed tolerances; incomplete integration is
%   an error unless a terminal Earth/outer guard event ended the trajectory.
%
% INPUT:
%   axis             [1 x 4] Magnitude [VU], azimuth/elevation [rad], phase [TU].
%   departure_state  [6 x 1] Unperturbed rotating barycentric state in DU/VU.
%   settings         Frozen model, units, event radii and integration settings;
%                    see data_builders/build_rpo_to_moon_database.m.
%
% OUTPUT:
%   rows     [H x 20] Historical encounter table; H may be zero.
%   retried  Logical flag indicating use of the relaxed-tolerance attempt.

% Apply exactly the historical rotating-frame Cartesian impulse.
dv = abs(axis(1));
theta = axis(2);
inclination = axis(3);
impulse = [dv*cos(theta)*cos(inclination), ...
    dv*sin(theta)*cos(inclination), dv*sin(inclination)];
initial_state = departure_state(:);
initial_state(4:6) = initial_state(4:6) + impulse.';
final_time = settings.direction * settings.horizon_days / settings.time_unit_days;
options = odeset('RelTol', settings.relative_tolerance, ...
    'AbsTol', settings.absolute_tolerance, 'Events', @(t, x) events(t, x, settings));
rhs = @(t, x) rhs_CR3BP_3D(t, x, settings.mu);
retried = false;
try
    solution = ode113(rhs, [0, final_time], initial_state, options);
    require_complete(solution, final_time);
catch first_error
    retried = true;
    options = odeset(options, 'RelTol', settings.relative_tolerance*settings.retry_factor, ...
        'AbsTol', settings.absolute_tolerance*settings.retry_factor);
    try
        solution = ode113(rhs, [0, final_time], initial_state, options);
        require_complete(solution, final_time);
    catch retry_error
        failure = MException('rpomoon:IntegrationFailure', ...
            'Departure [%.17g %.17g %.17g %.17g] failed both integrations.', ...
            axis(1), axis(2), axis(3), axis(4));
        failure = addCause(failure, first_error);
        failure = addCause(failure, retry_error);
        throw(failure);
    end
end

% Preserve column order/shapes even for a trajectory with no lunar entries.
hits = find(solution.ie == 1);
count = numel(hits);
rows = array2table(repmat(axis, count, 1), 'VariableNames', ...
    {'deltaV', 'thetaV', 'verticalThetaV', 'orbitPhase'});
rows.deltaV_vect = repmat(impulse, count, 1);
rows.v_CR3BP = zeros(count, 3);
rows.v_CR3BP_mag = zeros(count, 1);
rows.n_hits = (1:count).';
rows.t_event = zeros(count, 1);
rows.t_event_days = zeros(count, 1);
rows.v_sc_ECI = zeros(count, 3);
rows.vInfMoon_RTH = zeros(count, 3);
rows.vInfMagMoon = zeros(count, 1);
rows.thetaInfMoon = zeros(count, 1);
rows.phiInfMoon = zeros(count, 1);
rows.thetaMoon_deg = zeros(count, 1);
rows.r_RTH = zeros(count, 3);
rows.t_RTH = zeros(count, 3);
rows.h_RTH = zeros(count, 3);
rows.v_moon_ECI = zeros(count, 3);
velocity_unit = settings.distance_unit_km / (settings.time_unit_days*86400);

% Convert each event using the original planar ECI transform and RTH basis.
for index = 1:count
    time = solution.xe(hits(index));
    state = solution.ye(:, hits(index)).';
    velocity = state(4:6);
    spacecraft = rot2in(state([1, 2, 4, 5]), time, settings.mu);
    moon = rot2in([1-settings.mu, 0, 0, 0], time, settings.mu);
    spacecraft_velocity = [spacecraft(3:4), velocity(3)] * velocity_unit;
    moon_velocity = [moon(3:4), 0] * velocity_unit;
    moon_position = [moon(1:2), 0] * settings.distance_unit_km;
    radial = moon_position/norm(moon_position);
    normal = [0, 0, 1];
    transverse = cross(normal, radial);
    relative_velocity = spacecraft_velocity - moon_velocity;
    local_velocity = [dot(relative_velocity, radial), ...
        dot(relative_velocity, transverse), dot(relative_velocity, normal)];
    rows.v_CR3BP(index, :) = velocity;
    rows.v_CR3BP_mag(index) = norm(velocity);
    rows.t_event(index) = time;
    rows.t_event_days(index) = time * settings.time_unit_days;
    rows.v_sc_ECI(index, :) = spacecraft_velocity;
    rows.v_moon_ECI(index, :) = moon_velocity;
    rows.vInfMoon_RTH(index, :) = local_velocity;
    rows.vInfMagMoon(index) = norm(local_velocity);
    rows.thetaInfMoon(index) = atan2(local_velocity(2), local_velocity(1));
    rows.phiInfMoon(index) = atan2(local_velocity(3), hypot(local_velocity(1), local_velocity(2)));
    moon_phase = mod(time, 2*pi);
    if moon_phase == 0 && time > 0
        moon_phase = 2*pi; % Historical wrapTo2Pi positive-endpoint convention.
    end
    rows.thetaMoon_deg(index) = rad2deg(moon_phase);
    rows.r_RTH(index, :) = radial;
    rows.t_RTH(index, :) = transverse;
    rows.h_RTH(index, :) = normal;
end
rows.deltaV_vect(abs(rows.deltaV_vect) < 1e-10) = 0;
end

function [value, terminal, direction] = events(~, state, settings)
%EVENTS Lunar entry plus terminal Earth collision and outer escape guards.
mu = settings.mu;
earth_distance = sqrt((state(1)+mu)^2 + state(2)^2 + state(3)^2);
moon_distance = sqrt((state(1)-(1-mu))^2 + state(2)^2 + state(3)^2);
moon_boundary = settings.moon_radius_km/settings.distance_unit_km ...
    + settings.moon_soi_km/settings.distance_unit_km;
value = [moon_distance-moon_boundary; ...
    earth_distance-settings.earth_guard_km/settings.distance_unit_km; ...
    earth_distance-settings.outer_guard_km/settings.distance_unit_km];
terminal = [0; 1; 1];
direction = [-1; -1; 1];
end

function require_complete(solution, final_time)
%REQUIRE_COMPLETE Do not mistake an integrator's early warning for completion.
terminal_event = ~isempty(solution.ie) && any(solution.ie(end) == [2, 3]) ...
    && abs(solution.x(end)-solution.xe(end)) <= 1e-10*max(1, abs(final_time));
if any(~isfinite(solution.y(:))) || ...
        (solution.x(end) ~= final_time && ~terminal_event)
    error('rpomoon:IncompleteIntegration', ...
        'Integrator stopped at %.17g TU before the %.17g-TU horizon without a terminal guard.', ...
        solution.x(end), final_time);
end
end
