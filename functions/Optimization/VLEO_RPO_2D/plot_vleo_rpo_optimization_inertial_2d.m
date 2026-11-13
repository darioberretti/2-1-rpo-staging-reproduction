function plot_vleo_rpo_optimization_inertial_2d(setup, solution, figure_id, figure_name, axis_limits, title_suffix, show_title)
%PLOT_VLEO_RPO_OPTIMIZATION_INERTIAL_2D Draw a stored planar transfer in Earth-centered inertial coordinates.
%
% PROTOTYPE:
%   plot_vleo_rpo_optimization_inertial_2d(setup, solution, figure_id, figure_name, ...
%       axis_limits, title_suffix, show_title)
%
% DESCRIPTION:
%   Rotate the stored trajectory/nodes and align a two-period inertial RPO
%   reference with the arrival epoch. Positions remain normalized in DU.
%   Select and clear the supplied figure; the runner creates its docked export
%   container beforehand. Plotting does not execute a transfer optimization.
%
% INPUT:
%   setup       Multiple-shooting setup with constants, target orbit, and cfg.
%   solution    Postprocessed initial/optimized transfer with rotating DU/VU
%                 state rows and TU epochs.
%   figure_id   Existing figure number selected and cleared for drawing.
%   figure_name Descriptive name used for the figure/title.
%   axis_limits Optional [xmin, xmax, ymin, ymax] in DU; default [] is automatic.
%   title_suffix  Optional extra title text, default empty.
%   show_title    Optional logical, default true.
%
% OUTPUT:
%   None. Draws into figure_id; saving/numbering belongs to the runner exporter.

if nargin < 5
    axis_limits = [];
end

% Use no extra title suffix when the caller omits it.
if nargin < 6
    title_suffix = '';
end

% Show plot titles by default when the caller omits the flag.
if nargin < 7
    show_title = true;
end

% Convert the optimized trajectory into the Earth-centered inertial frame.
X_traj_in = rot2in(solution.X_traj, solution.t_traj, setup.constants.mu);

% Convert the optimized nodes into the Earth-centered inertial frame.
X_nodes_in = rot2in(solution.X_nodes, solution.dec.t_nodes, setup.constants.mu);

% Build one rotating-frame 3BP RPO reference orbit over one period only.
[~, ~, ~, x_rpo_rot, t_rpo] = propagate(setup.x0_rpo_planar_3bp, 0, setup.constants.T_RPO_3bp, setup.constants.mu, '2D', 'CR3BP', 'false', false, [], [], setup.cfg.abs_tol, setup.cfg.rel_tol);

% Shift the inertial RPO epoch so the point at tau = phase_rpo matches the transfer arrival epoch.
t_rpo_shift = solution.dec.t_nodes(end) - solution.phase_rpo;

% Rotate the same one-period orbit with the aligned and one-period-shifted inertial times.
X_rpo_in_first = rot2in(x_rpo_rot, t_rpo + t_rpo_shift, setup.constants.mu);
X_rpo_in_second = rot2in(x_rpo_rot, t_rpo + t_rpo_shift + setup.constants.T_RPO_3bp, setup.constants.mu);
X_rpo_in = [X_rpo_in_first; X_rpo_in_second(2:end, :)];

% Build one inertial departure parking circle at the initial epoch.
theta_circle = linspace(0, 2 * pi, 400);
r_depart_du = (setup.constants.R_e + solution.departure_altitude_km) / setup.constants.DU;
x_circle = r_depart_du * cos(theta_circle);
y_circle = r_depart_du * sin(theta_circle);

% Build one compact title with the main solution metrics.
title_text = sprintf('%s inertial%s | h = %.1f km | total DV = %.1f m/s', figure_name, title_suffix, solution.departure_altitude_km, solution.dv_total_mps);

% Create the inertial figure and plot all main geometry.
figure(figure_id);
clf;
set(gcf, 'Name', [figure_name ' inertial'], 'Theme', 'light');
hold on;
grid on;
box on;
plot(X_rpo_in(:, 1), X_rpo_in(:, 2), '-', 'Color', [0.3 0.5 0.9], 'LineWidth', 1.0, 'DisplayName', 'RPO');
plot(x_circle, y_circle, '-', 'Color', [0.85 0.2 0.2], 'LineWidth', 1.0, 'DisplayName', 'LEO circle', 'HandleVisibility', 'off');
plot(0, 0, '.', 'Color', [0.0 0.2 0.8], 'MarkerSize', 20, 'DisplayName', 'Earth');
plot(X_traj_in(:, 1), X_traj_in(:, 2), 'k', 'LineWidth', 1, 'DisplayName', 'Trajectory');
plot(X_nodes_in(:, 1), X_nodes_in(:, 2), 'o', 'Color', 'k', 'MarkerFaceColor', 'none', 'MarkerSize', 7, 'LineWidth', 0.7, 'DisplayName', 'Nodes');

% Label and format the inertial axes.
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
