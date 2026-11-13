function figures = plot_method_exploration(result, visible)
%PLOT_METHOD_EXPLORATION Plot the four Run 02 method diagnostics.
%
% PROTOTYPE:
%   figures = rpo2ds.plot_method_exploration(result, visible)
%
% DESCRIPTION:
%   Show the encounter map, Earth-departure construction, heliocentric
%   target transfers, and perigee-altitude comparison using stored data.
%
% INPUT:
%   result   Study result supplying data.encounter_map, earth_departure, target_transfers,
%              altitude_sensitivity, targets, and configuration.shared.
%   visible  Figure visibility, on or off.
%
% OUTPUT:
%   figures  Struct of registered names and figure handles, ready for export.
%              Plotting does not save result files or execute optimizations.

figures.rpo_moon_encounter_map = local_encounter_map(result, visible);
figures.earth_departure_construction = local_earth_departure(result, visible);
figures.venus_mars_target_transfers = local_target_transfers(result, visible);
figures.escape_direction_altitude_comparison = ...
    local_direction_comparison(result, visible);
end

function handle = local_encounter_map(result, visible)
%LOCAL_ENCOUNTER_MAP Plot the deterministic encounter subset colored by Moon-SOI hit count.

data = result.data.encounter_map;
rows = data.source_rows;
handle = rpo.create_figure(visible,[900,620]);
scatter(data.time_of_flight_days(rows), data.vinf_moon_km_s(rows), 10, ...
    data.moon_soi_entries(rows), 'filled');
grid on; box on;
xlabel('RPO departure-to-Moon encounter TOF [d]', 'Interpreter', 'latex');
ylabel('$|v_{\infty,M}^{-}|$ [km/s]', 'Interpreter', 'latex');
cb = colorbar;
cb.Label.String = 'Moon SOI entries';
cb.Ticks = unique(data.moon_soi_entries(rows), 'sorted');
colormap(handle, parula(256));
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
end

function handle = local_earth_departure(result, visible)
%LOCAL_EARTH_DEPARTURE Draw perigee geometry and representative Moon-Earth Lambert arcs.

data = result.data.earth_departure;
shared = result.configuration.shared;
handle = rpo.create_figure(visible,[1150,570]);
tiledlayout(handle, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on; grid on; box on; axis equal; view(45, 20);
[sx, sy, sz] = sphere(30);
earth_radius = shared.earth_radius_km;
surf(earth_radius * sx, earth_radius * sy, earth_radius * sz, ...
    'FaceColor', [0.2, 0.45, 0.9], 'EdgeColor', 'none', ...
    'FaceAlpha', 0.85, 'DisplayName', 'Earth');
indices = round(linspace(1, size(data.perigee_locus_km, 1), 21));
for index = 1:numel(indices)
    row = indices(index);
    arc = local_two_body_path(data.perigee_locus_km(row, :), ...
        data.perigee_velocity_km_s(row, :), linspace(0, 6000, 120), ...
        shared.earth_gravitational_parameter_km3_s2);
    if index == 1
        plot3(arc(:, 1), arc(:, 2), arc(:, 3), ...
            'Color', [0.8, 0, 0.7], 'LineWidth', 0.8, ...
            'DisplayName', 'Outgoing hyperbolae');
    else
        plot3(arc(:, 1), arc(:, 2), arc(:, 3), ...
            'Color', [0.8, 0, 0.7], 'LineWidth', 0.8, ...
            'HandleVisibility', 'off');
    end
end
plot3(data.perigee_locus_km(:, 1), data.perigee_locus_km(:, 2), ...
    data.perigee_locus_km(:, 3), 'r-', 'LineWidth', 2, ...
    'DisplayName', 'Earth-perigee locus');
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
title('Admissible Earth-perigee locus');
legend('Location', 'northeast');
xlim([-30000, 30000]); ylim([-30000, 30000]); zlim([-30000, 30000]);

nexttile; hold on; grid on; box on; axis equal; view(-45, 25);
moon_radius = shared.distance_unit_km;
angle = linspace(0, 2 * pi, 300);
plot3(moon_radius * cos(angle), moon_radius * sin(angle), ...
    zeros(size(angle)), '--', 'Color', [0.55, 0.55, 0.55], ...
    'DisplayName', 'Moon orbit');
colors = turbo(numel(data.lambert_arcs));
for index = 1:numel(data.lambert_arcs)
    arc = data.lambert_arcs{index};
    plot3(arc.positions_km(:, 1), arc.positions_km(:, 2), arc.positions_km(:, 3), ...
        'Color', colors(index, :), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%.1f d', arc.time_of_flight_days));
    moon_marker = plot3(arc.moon_position_km(1), arc.moon_position_km(2), ...
        arc.moon_position_km(3), 'ko', 'MarkerSize', 5, ...
        'MarkerFaceColor', 'k', 'DisplayName', 'Moon departure');
    if index > 1, moon_marker.HandleVisibility = 'off'; end
end
plot3(0, 0, 0, 'ko', 'MarkerFaceColor', [0.2, 0.5, 0.9], 'DisplayName', 'Earth');
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
title('Moon-to-Earth Lambert arcs');
legend('Location', 'northeast');
set(findall(handle, '-property', 'FontName'), 'FontName', 'Times New Roman');
end

function handle = local_target_transfers(result, visible)
%LOCAL_TARGET_TRANSFERS Draw the stored heliocentric transfers and planetary orbit context.

% The tracked compact target input supplies the accepted planetary states.
repository_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
targets = rpo2ds.load_targets(repository_root);
mu_sun = targets.sun_gravitational_parameter_km3_s2;
handle = rpo.create_figure(visible,[1150,560]);
tiledlayout(handle, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
names = ["venus", "mars"];
display_names = ["Venus", "Mars"];
for index = 1:2
    transfer = result.data.target_transfers.(names(index));
    target = targets.(names(index));
    earth_orbit = local_two_body_orbit(target.earth_departure_position_km, ...
        target.earth_departure_velocity_km_s, mu_sun);
    planet_orbit = local_two_body_orbit(target.planet_arrival_position_km, ...
        target.planet_arrival_velocity_km_s, mu_sun);
    ax = nexttile;
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on'); axis(ax, 'equal');
    view(ax, -45, 25);
    plot3(earth_orbit(:, 1), earth_orbit(:, 2), earth_orbit(:, 3), ...
        'b-', 'LineWidth', 0.9, 'DisplayName', 'Earth orbit');
    plot3(planet_orbit(:, 1), planet_orbit(:, 2), planet_orbit(:, 3), ...
        'r-', 'LineWidth', 0.9, ...
        'DisplayName', char(display_names(index) + " orbit"));
    plot3(transfer.positions_km(:, 1), transfer.positions_km(:, 2), ...
        transfer.positions_km(:, 3), 'k-', 'LineWidth', 1.3, ...
        'DisplayName', 'Transfer');
    plot3(transfer.earth_position_km(1), transfer.earth_position_km(2), ...
        transfer.earth_position_km(3), 'bo', 'MarkerFaceColor', 'b', ...
        'HandleVisibility', 'off');
    plot3(transfer.planet_position_km(1), transfer.planet_position_km(2), ...
        transfer.planet_position_km(3), 'ro', 'MarkerFaceColor', 'r', ...
        'HandleVisibility', 'off');
    plot3(0, 0, 0, 'o', 'Color', [0.9, 0.6, 0], ...
        'MarkerFaceColor', [0.9, 0.6, 0], 'DisplayName', 'Sun');
    limit = 1.12 * max([vecnorm(earth_orbit, 2, 2); vecnorm(planet_orbit, 2, 2)]);
    xlim(ax, limit * [-1, 1]); ylim(ax, limit * [-1, 1]);
    zlim(ax, 0.5 * limit * [-1, 1]);
    xlabel(ax, 'x [km]'); ylabel(ax, 'y [km]'); zlabel(ax, 'z [km]');
    title(ax, display_names(index));
    legend(ax, 'Location', 'northeast');
    set(ax, 'FontName', 'Times New Roman', 'FontSize', 13);
end
end

function positions = local_two_body_orbit(position, velocity, mu)
%LOCAL_TWO_BODY_ORBIT Propagate one two-body orbital period for the context curve.

semimajor_axis = 1 / (2 / norm(position) - dot(velocity, velocity) / mu);
period = 2 * pi * sqrt(semimajor_axis^3 / mu);
times = linspace(0, period, 300);
positions = local_two_body_path(position, velocity, times, mu);
end

function positions = local_two_body_path(position, velocity, times, mu)
%LOCAL_TWO_BODY_PATH Return inertial position rows at the supplied output times.

[~, states] = ode113(@(~, state) [state(4:6); ...
    -mu * state(1:3) / norm(state(1:3))^3], times, ...
    [position(:); velocity(:)], odeset('RelTol', 1e-10, 'AbsTol', 1e-11));
positions = states(:, 1:3);
end

function handle = local_direction_comparison(result, visible)
%LOCAL_DIRECTION_COMPARISON Plot the two escape-direction costs against Earth-perigee altitude.

comparison = result.data.altitude_sensitivity.direction_comparison;
altitudes = result.data.altitude_sensitivity.altitude_km;
labels = {'In-plane (x+y)', 'Out of plane (+z)'};
colors = [0.15, 0.35, 0.8; 0.55, 0.15, 0.55];
handle = rpo.create_figure(visible,[1050,720]);
tiledlayout(handle, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile; hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');
for index = 1:numel(labels)
    plot(ax1, altitudes, comparison.best_delta_v_km_s(:, index), '-o', ...
        'Color', colors(index, :), 'LineWidth', 1.7, 'MarkerSize', 5, ...
        'DisplayName', labels{index});
end
set(ax1, 'XScale', 'log');
xlim(ax1, [min(altitudes), max(altitudes)]);
ylabel(ax1, 'Best patched Earth-perigee correction [km/s]');
title(ax1, sprintf('Fixed Earth-escape speed: %.0f km/s', comparison.speed_km_s));
legend(ax1, 'Location', 'northwest', 'AutoUpdate', 'off');
xline(ax1, 1000, '--k');
first_vertical_seed = altitudes(find(isfinite( ...
    comparison.best_delta_v_km_s(:, 2)), 1));
text(ax1, 1300, 0.82 * max(comparison.best_delta_v_km_s(:), [], 'omitnan'), ...
    sprintf('First sampled +z seed: %.0f km', first_vertical_seed), ...
    'Color', colors(2, :), 'FontSize', 11);

ax2 = nexttile; hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');
for index = 1:numel(labels)
    plot(ax2, altitudes, ...
        comparison.n_candidates_in_rpo_vinf_range(:, index), '-o', ...
        'Color', colors(index, :), 'LineWidth', 1.7, 'MarkerSize', 5);
end
set(ax2, 'XScale', 'log');
xlim(ax2, [min(altitudes), max(altitudes)]);
xline(ax2, 1000, '--k');
text(ax2, 1200, 0.9 * max(comparison.n_candidates_in_rpo_vinf_range(:)), ...
    'Search starts at 1000 km', 'FontSize', 11);
xlabel(ax2, 'Earth-perigee altitude [km]');
ylabel(ax2, 'Candidates after lunar-v_\infty range filter');
sgtitle(handle, 'Why the patched-seed search varies Earth-perigee altitude', ...
    'FontName', 'Times New Roman');
set([ax1, ax2], 'FontName', 'Times New Roman', 'FontSize', 13);
end
