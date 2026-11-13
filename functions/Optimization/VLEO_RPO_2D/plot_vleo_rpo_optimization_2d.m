function plot_vleo_rpo_optimization_2d(setup, solution, figure_id, figure_name, show_saved_seed_reference, axis_limits, title_suffix, show_title)
%PLOT_VLEO_RPO_OPTIMIZATION_2D Draw a stored planar transfer in rotating barycentric coordinates.
%
% PROTOTYPE:
%   plot_vleo_rpo_optimization_2d(setup, solution, figure_id, figure_name, ...
%       show_saved_seed_reference, axis_limits, title_suffix, show_title)
%
% DESCRIPTION:
%   Overlay the saved seed when requested, a nominal CR3BP reference RPO,
%   departure circle, active impulses, trajectory, and shooting nodes in DU.
%   Select and clear the supplied figure; the runner creates its docked export
%   container beforehand. Plotting does not execute a transfer optimization.
%
% INPUT:
%   setup       Multiple-shooting setup with constants, target orbit, and cfg.
%   solution    Postprocessed initial/optimized transfer with rotating DU/VU
%                 state rows and TU epochs.
%   figure_id   Existing figure number selected and cleared for drawing.
%   figure_name Descriptive name used for the figure/title.
%   show_saved_seed_reference  Optional logical, default true.
%   axis_limits Optional [xmin, xmax, ymin, ymax] in DU; default [] is automatic.
%   title_suffix  Optional extra title text, default empty.
%   show_title    Optional logical, default true.
%
% OUTPUT:
%   None. Draws into figure_id; saving/numbering belongs to the runner exporter.

if nargin < 5
    show_saved_seed_reference = true;
end

% Use automatic axes when the caller omits the close-up limits.
if nargin < 6
    axis_limits = [];
end

% Use no extra title suffix when the caller omits it.
if nargin < 7
    title_suffix = '';
end

% Show plot titles by default when the caller omits the flag.
if nargin < 8
    show_title = true;
end

% Build the figure that hosts the trajectory plot.
figure(figure_id);
clf;
set(gcf, 'Name', figure_name, 'Theme', 'light');
hold on;
grid on;
box on;

% Read the main constants used by the plot.
mu = setup.constants.mu;
DU = setup.constants.DU;
title_text = sprintf('%s rot%s | h = %.1f km | total DV = %.1f m/s', figure_name, title_suffix, solution.departure_altitude_km, solution.dv_total_mps);

% Plot the saved continuous seed only when requested.
if show_saved_seed_reference
    plot(setup.seed_state_rot(:, 1), setup.seed_state_rot(:, 2), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.0, 'DisplayName', 'Saved seed');
end

% Plot the nominal one-period 3BP RPO target.
[~, ~, ~, x_rpo_plot, ~] = propagate(setup.x0_rpo_planar_3bp, 0, setup.constants.T_RPO_3bp, setup.constants.mu, '2D', 'CR3BP', 'false', false, [], [], setup.cfg.abs_tol, setup.cfg.rel_tol);
plot(x_rpo_plot(:, 1), x_rpo_plot(:, 2), '-', 'Color', [0.3 0.5 0.9], 'LineWidth', 1.0, 'DisplayName', 'RPO');

% Plot the departure parking circle at the active departure altitude.
theta_circle = linspace(0, 2 * pi, 400);
r_depart_du = (setup.constants.R_e + solution.departure_altitude_km) / DU;
plot(-mu + r_depart_du * cos(theta_circle), r_depart_du * sin(theta_circle), '-', 'Color', [0.85 0.2 0.2], 'LineWidth', 1.0, 'DisplayName', 'LEO circle', 'HandleVisibility', 'off');

% Plot the Earth and the Moon in the rotating frame.
plot(-mu, 0, '.', 'Color', [0.0 0.2 0.8], 'MarkerSize', 20, 'DisplayName', 'Earth');
plot(1 - mu, 0, '.', 'Color', [0.4 0.4 0.4], 'MarkerSize', 18, 'DisplayName', 'Moon');

% Plot the active node impulses as quivers.
active_rows = find(setup.active_impulse_mask);
quiver(solution.X_nodes(active_rows, 1), solution.X_nodes(active_rows, 2), solution.DV_nodes_rot(active_rows, 1), solution.DV_nodes_rot(active_rows, 2), 0, 'Color', [0.9 0.2 0.2], 'LineWidth', 1.0, 'DisplayName', 'Impulses');

% Plot the current continuous trajectory last so it stays above the RPO curve.
plot(solution.X_traj(:, 1), solution.X_traj(:, 2), 'k', 'LineWidth', 1, 'DisplayName', 'Trajectory');

% Plot the current shooting nodes after the trajectory so they remain visible.
plot(solution.X_nodes(:, 1), solution.X_nodes(:, 2), 'o', 'Color', 'k', 'MarkerFaceColor', 'none', 'MarkerSize', 7, 'LineWidth', 0.7, 'DisplayName', 'Nodes');

% Label the axes and optionally draw the plot title.
xlabel('x [DU]');
ylabel('y [DU]');
if show_title
    title(title_text, 'Interpreter', 'none');
end
if isempty(axis_limits)
    axis equal;
else
    axis equal;
    xlim(axis_limits(1:2));
    ylim(axis_limits(3:4));
end
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Location', 'best');

end
