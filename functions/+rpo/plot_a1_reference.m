function figures = plot_a1_reference(result, visible)
%PLOT_A1_REFERENCE Plot the A1 samples and corrected reference RPO.
%
% PROTOTYPE:
%   figures = rpo.plot_a1_reference(result, visible)
%
% DESCRIPTION:
%   Create the sampled-family overview and the reference-orbit overlay.
%
% INPUT:
%   result   Finalized study result supplying
%              data.figure_data, data.rpo.trajectory, and configuration.constants.
%   visible  Figure visibility, on or off.
%
% OUTPUT:
%   figures  Struct of registered figure names and handles, ready for
%              rpo.export_figures. Plotting does not write files.

mu = result.configuration.constants.mu_earth_moon;
samples = result.data.figure_data.family_sample_trajectories;
sample_indices = result.data.figure_data.family_sample_indices;
sample_rows = result.data.figure_data.family_sample_rows;

% Show the frozen family samples in individual panels.
figures.a1_family_samples = rpo.create_figure(visible,[1200,800]);
columns = min(5, numel(samples));
rows = ceil(numel(samples) / columns);
tiledlayout(rows, columns, 'TileSpacing', 'compact', 'Padding', 'compact');
for index = 1:numel(samples)
    nexttile;
    trajectory = samples{index};
    plot(trajectory(:, 1), trajectory(:, 2), 'k', 'LineWidth', 1.1);
    hold on;
    plot(-mu, 0, 'k.', 'MarkerSize', 9);
    plot(1 - mu, 0, 'k.', 'MarkerSize', 9);
    title(sprintf('$n = %d, x_0 = %.2f$', ...
        sample_indices(index), sample_rows(index, 1)), ...
        'Interpreter', 'latex', 'FontSize', 10);
    xlabel('$x$ [DU]', 'Interpreter', 'latex');
    ylabel('$y$ [DU]', 'Interpreter', 'latex');
    axis equal;
    grid on;
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 9);
end

% Overlay the corrected reference orbit on the selected family trajectories.
figures.reference_rpo_family = rpo.create_figure(visible,[]);
hold on;
for index = 1:result.data.figure_data.rpo_neighbor_count
    trajectory = samples{index};
    plot(trajectory(:, 1), trajectory(:, 2), 'k', 'LineWidth', 1, ...
        'HandleVisibility', 'off');
end
trajectory = result.data.rpo.trajectory;
plot(trajectory(:, 1), trajectory(:, 2), 'r', 'LineWidth', 2.5, ...
    'DisplayName', '2:1 RPO');
plot(-mu, 0, 'k.', 'MarkerSize', 18, 'HandleVisibility', 'off');
plot(1 - mu, 0, 'k.', 'MarkerSize', 18, 'HandleVisibility', 'off');
text(-mu, -0.10, 'Earth', 'Interpreter', 'latex', ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(1.18, -0.05, 'Moon', 'Interpreter', 'latex', ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold');
xlabel('$x$ [DU]', 'Interpreter', 'latex');
ylabel('$y$ [DU]', 'Interpreter', 'latex');
title('2:1 synodic-resonance RPO', 'Interpreter', 'latex');
legend('show', 'Interpreter', 'latex', 'Location', 'best');
axis equal;
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
end
