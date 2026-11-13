function figures = plot_validation(result, target_name, visible)
%PLOT_VALIDATION Create the four validation figures for one target.
%
% PROTOTYPE:
%   figures = rpo2ds.plot_validation(result, target_name, visible)
%
% DESCRIPTION:
%   Show the patched initial guess, shooting nodes, and optimized transfer
%   in rotating DU and inertial km. Impulse arrows are scaled for display.
%
% INPUT:
%   result   Study result supplying data from pack_transfer_result and configuration.shared.
%   target_name  venus or mars; selects the registered figure-name prefix.
%   visible  Figure visibility, on or off.
%
% OUTPUT:
%   figures  Struct of registered names and figure handles, ready for export.
%              Plotting does not save result files or execute optimizations.

figures = struct();
figures.([target_name '_patched_initial_guess']) = ...
    local_initial_guess(result.data, result.configuration.shared, visible);
figures.([target_name '_multiple_shooting_nodes']) = ...
    local_nodes(result.data, visible);
figures.([target_name '_optimized_rotating']) = ...
    local_optimized_rotating(result.data, visible);
figures.([target_name '_optimized_inertial']) = ...
    local_optimized_inertial(result.data, result.configuration.shared, visible);
end

function handle = local_initial_guess(data, shared, visible)
%LOCAL_INITIAL_GUESS Plot the patched rotating-frame arcs, Moon SOI, and terminal points.

handle = local_figure(visible, [900, 650]);
hold on; grid on; box on; axis equal;
initial = data.initial_guess;
mu = shared.earth_moon_mass_parameter;
soi_radius = (mu / (1 - mu))^(2/5);
angle = linspace(0, 2 * pi, 200);
fill(1 - mu + soi_radius * cos(angle), soi_radius * sin(angle), ...
    [1, 1, 0.6], 'EdgeColor', 'none', 'DisplayName', 'Moon SOI');
plot(initial.xx(:, 1), initial.xx(:, 2), 'Color', [0, 0.447, 0.741], ...
    'LineWidth', 1.3, 'DisplayName', 'RPO--Moon leg');
plot(initial.xxRot(:, 1), initial.xxRot(:, 2), ...
    'Color', [0.85, 0.325, 0.098], 'LineWidth', 1.3, ...
    'DisplayName', 'Moon--Earth leg');
plot(initial.xxFB(:, 1), initial.xxFB(:, 2), ...
    'Color', [0.929, 0.694, 0.125], 'LineWidth', 1.3, ...
    'DisplayName', 'Moon flyby leg');
plot(initial.xx(1, 1), initial.xx(1, 2), 'bo', ...
    'MarkerSize', 10, 'LineWidth', 2.5, 'DisplayName', 'Initial point');
plot(initial.xxRot(end, 1), initial.xxRot(end, 2), 'ro', ...
    'MarkerSize', 10, 'LineWidth', 2.5, 'DisplayName', 'Final point');
xlabel('x [DU]'); ylabel('y [DU]');
legend('Location', 'northeastoutside');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
end

function handle = local_nodes(data, visible)
%LOCAL_NODES Overlay available initial shooting nodes and their regularized impulse seed.

handle = local_figure(visible, [900, 650]);
hold on; grid on; box on; axis equal;
nodes = data.nodes;
trajectory = local_positions(data.initial_guess, 'xxTot');
if ~isempty(trajectory)
    plot(trajectory(:, 1), trajectory(:, 2), 'k-', ...
        'LineWidth', 1.1, 'HandleVisibility', 'off');
end
node_points = local_positions(nodes, 'X_nodes_IG');
if isempty(node_points)
    node_points = local_positions(nodes, 'xx_nodes_IG');
end
if ~isempty(node_points)
    plot(node_points(:, 1), node_points(:, 2), 'ko', ...
        'MarkerSize', 8, 'DisplayName', 'Nodes');
end
if isfield(data.initial_guess, 'nlp')
    guess = data.initial_guess.nlp;
    local_impulses(guess.xx_nodes, guess.DVvec_opt_rot, trajectory);
end
if isempty(trajectory) && isempty(node_points)
    local_placeholder('Node-layout data unavailable');
else
    xlabel('x [DU]'); ylabel('y [DU]');
    legend('Location', 'northeastoutside', 'Interpreter', 'latex');
end
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
end

function handle = local_optimized_rotating(data, visible)
%LOCAL_OPTIMIZED_ROTATING Plot stored optimized rotating-frame positions and scaled burn arrows.

handle = local_figure(visible, [900, 650]);
hold on; grid on; box on; axis equal;
trajectory = local_positions(data.optimized, 'xx_traj_rot');
nodes = local_positions(data.optimized, 'xx_nodes');
if ~isempty(trajectory)
    plot(trajectory(:, 1), trajectory(:, 2), 'k-', ...
        'LineWidth', 1.2, 'HandleVisibility', 'off');
end
if ~isempty(nodes)
    plot(nodes(:, 1), nodes(:, 2), 'ko', ...
        'MarkerSize', 8, 'DisplayName', 'Nodes');
end
local_impulses(nodes, data.optimized_impulses.rotating_nd, trajectory);
if isempty(trajectory)
    local_placeholder('Optimized rotating-frame data unavailable');
else
    xlabel('x [DU]'); ylabel('y [DU]');
    legend('Location', 'northeastoutside', 'Interpreter', 'latex');
end
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
end

function handle = local_optimized_inertial(data, shared, visible)
%LOCAL_OPTIMIZED_INERTIAL Plot stored inertial positions [km], the Moon path, and scaled burns.

handle = local_figure(visible, [900, 650]);
hold on; grid on; box on; axis equal;
trajectory = local_positions(data.optimized, 'xx_traj_in');
nodes = local_positions(data.optimized, 'xx_nodes_in');
if ~isempty(trajectory)
    plot(trajectory(:, 1), trajectory(:, 2), 'k-', ...
        'LineWidth', 1.2, 'HandleVisibility', 'off');
end
if ~isempty(nodes)
    plot(nodes(:, 1), nodes(:, 2), 'ko', ...
        'MarkerSize', 8, 'DisplayName', 'Nodes');
end
angle = linspace(0, 2 * pi, 200);
moon_radius = shared.distance_unit_km;
plot(moon_radius * cos(angle), moon_radius * sin(angle), '--', ...
    'Color', [0.6, 0.6, 0.6], 'DisplayName', 'Moon path');
local_impulses(nodes, data.optimized_impulses.inertial_m_s, trajectory);
if isempty(trajectory)
    local_placeholder('Optimized inertial-frame data unavailable');
else
    xlabel('x [km]'); ylabel('y [km]');
    legend('Location', 'northeastoutside', 'Interpreter', 'latex');
end
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
end

function local_impulses(nodes, vectors, trajectory)
%LOCAL_IMPULSES Scale burn arrows to the displayed trajectory span and suppress tiny burns.

if isempty(nodes) || isempty(vectors) || isempty(trajectory)
    return
end
strength = vecnorm(vectors(:, 1:2), 2, 2);
active = strength > 1e-4 * max(strength);
span = max(range(trajectory(:, 1)), range(trajectory(:, 2)));
scale = 0.38 * span / max(strength);
quiver(nodes(active, 1), nodes(active, 2), ...
    scale * vectors(active, 1), scale * vectors(active, 2), 0, ...
    'Color', 'b', 'LineWidth', 0.8, 'DisplayName', '$\Delta V$ impulses');
end

function handle = local_figure(visible, position)
%LOCAL_FIGURE Create a docked figure with the requested export width/height [pixels].

handle = rpo.create_figure(visible,position);
end

function points = local_positions(source, field_name)
%LOCAL_POSITIONS Extract the first three numeric columns, or [] when position data are absent.

points = [];
if ~isfield(source, field_name) || isempty(source.(field_name))
    return
end
value = source.(field_name);
if ~isnumeric(value) || size(value, 2) < 3
    return
end
points = value(:, 1:3);
end

function local_placeholder(message)
%LOCAL_PLACEHOLDER Place a missing-data message at the center of the current axes.

text(0.5, 0.5, message, 'Units', 'normalized', 'HorizontalAlignment', 'center');
end
