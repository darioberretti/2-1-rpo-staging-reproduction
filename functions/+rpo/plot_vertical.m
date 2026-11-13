function figures = plot_vertical(result, visible)
%PLOT_VERTICAL Plot vertical-magnitude, trajectory, and phase studies.
%
% PROTOTYPE:
%   figures = rpo.plot_vertical(result, visible)
%
% DESCRIPTION:
%   Create signed-magnitude time histories, one spatial trajectory, and
%   phase-dependent time histories. Times are in days and positions in DU.
%
% INPUT:
%   result   Finalized study result supplying
%              data.magnitude, data.phase, and configuration.
%   visible  Figure visibility, on or off.
%
% OUTPUT:
%   figures  Struct of registered figure names and handles, ready for
%              rpo.export_figures. Plotting does not write files.

impulses = result.data.magnitude.impulses_vu;
magnitudes = unique(abs(impulses), 'stable');

% Compare positive and negative impulses at each selected magnitude.
figures.rpo_vertical_perturbation_magnitudes = rpo.create_figure(visible,[900,900]);
tiledlayout(numel(magnitudes), 1, 'TileSpacing', 'compact');
for magnitude_index = 1:numel(magnitudes)
    nexttile;
    hold on;
    selected_impulses = find(abs(impulses) == magnitudes(magnitude_index));
    for impulse = selected_impulses'
        time = result.data.magnitude.times_days{impulse};
        trajectory = result.data.magnitude.trajectories{impulse};
        plot(time, trajectory(:, 3), 'LineWidth', 1.4);
    end
    ylabel('$z$ [DU]', 'Interpreter', 'latex');
    title(sprintf('$|v_{z0}| = %s$ VU', ...
        local_impulse_label(magnitudes(magnitude_index))), ...
        'Interpreter', 'latex');
    ylim([-1, 1]);
    grid on;
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
end
xlabel('$t$ [days]', 'Interpreter', 'latex');
legend({'Negative impulse', 'Positive impulse'}, 'Interpreter', 'latex');

% Show the first phase-case trajectory with its start/end points and primaries.
phase_impulse = result.data.phase.impulse_vu;
trajectory = result.data.phase.trajectories{1};
figures.rpo_vertical_trajectory = rpo.create_figure(visible,[]);
plot3(trajectory(:, 1), trajectory(:, 2), trajectory(:, 3), ...
    'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.2);
hold on;
plot3(trajectory(1, 1), trajectory(1, 2), trajectory(1, 3), ...
    'r.', 'MarkerSize', 18);
plot3(trajectory(end, 1), trajectory(end, 2), trajectory(end, 3), ...
    'r.', 'MarkerSize', 18);
mu = result.configuration.constants.mu_earth_moon;
plot3(-mu, 0, 0, 'k.', 'MarkerSize', 14);
plot3(1 - mu, 0, 0, 'k.', 'MarkerSize', 14);
text(trajectory(1, 1), trajectory(1, 2), trajectory(1, 3), ...
    ' Start', 'Color', 'r');
text(trajectory(end, 1), trajectory(end, 2), trajectory(end, 3), ...
    ' End', 'Color', 'r');
text(-mu, 0, 0, ' Earth');
text(1 - mu, 0, 0, ' Moon');
xlabel('$x$ [DU]', 'Interpreter', 'latex');
ylabel('$y$ [DU]', 'Interpreter', 'latex');
zlabel('$z$ [DU]', 'Interpreter', 'latex');
title(sprintf('%.0f days propagation ($v_{z0}=%s$ VU)', ...
    result.configuration.reference_rpo.vertical.duration_days, ...
    local_impulse_label(phase_impulse)), ...
    'Interpreter', 'latex');
axis equal;
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);

% Color the vertical-response histories by the initial orbital phase fraction.
figures.rpo_vertical_phase_response = rpo.create_figure(visible,[]);
hold on;
colors = copper(numel(result.data.phase.fractions));
for phase = 1:numel(result.data.phase.fractions)
    time = result.data.phase.propagation_times_days{phase};
    trajectory = result.data.phase.trajectories{phase};
    plot(time, trajectory(:, 3), 'Color', colors(phase, :), 'LineWidth', 1.3);
end
colormap(colors);
phase_limits = [min(result.data.phase.fractions), ...
    max(result.data.phase.fractions)];
if phase_limits(1) == phase_limits(2)
    phase_limits(2) = phase_limits(1) + 1;
end
clim(phase_limits);
colorbar_handle = colorbar;
colorbar_handle.Ticks = result.data.phase.fractions;
colorbar_handle.TickLabels = compose('%.1fT', result.data.phase.fractions);
xlabel('$t$ [days]', 'Interpreter', 'latex');
ylabel('$z$ [DU]', 'Interpreter', 'latex');
title(sprintf('$z$ response to $v_{z0}=%s$ VU at %d RPO phases', ...
    local_impulse_label(phase_impulse), numel(result.data.phase.fractions)), ...
    'Interpreter', 'latex');
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
end

function label = local_impulse_label(value)
%LOCAL_IMPULSE_LABEL Format an impulse magnitude for a LaTeX plot label.
%
% Use a power of ten for exact decades; otherwise use two significant digits.

exponent = round(log10(abs(value)));
if abs(abs(value) - 10 ^ exponent) <= 10 * eps(abs(value))
    label = sprintf('10^{%d}', exponent);
else
    label = sprintf('%.2g', value);
end
end
