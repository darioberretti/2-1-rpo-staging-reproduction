function figures = plot_reachability(result, visible)
%PLOT_REACHABILITY Plot the directional curve and symmetry reconstruction.
%
% PROTOTYPE:
%   figures = rpo2ds.plot_reachability(result, visible)
%
% DESCRIPTION:
%   Interpolate the converged polar sweep and reconstruct its historical
%   symmetric speed surface for display. Surface samples are not additional
%   optimized directions.
%
% INPUT:
%   result   Study result supplying data polar angles [rad], speeds [km/s], converged flags,
%              and renderer sampling settings.
%   visible  Figure visibility, on or off.
%
% OUTPUT:
%   figures  Struct of registered names and figure handles, ready for export.
%              Plotting does not save result files or execute optimizations.

figures.escape_reachability_curve = local_curve(result, visible);
figures.escape_velocity_cloud = local_cloud(result, visible);
end

function handle = local_curve(result, visible)
%LOCAL_CURVE Plot converged polar samples and their interpolated speed curve.

data = result.data;
handle = rpo.create_figure(visible,[850,560]);
valid = data.converged & isfinite(data.v_inf_magnitude_km_s);
phi = data.polar_angle_rad(valid);
speed = data.v_inf_magnitude_km_s(valid);
[phi, order] = sort(phi);
speed = speed(order);
fine_phi = linspace(phi(1), phi(end), data.renderer.polar_samples);
fine_speed = interp1(phi, speed, fine_phi, 'pchip');
plot(fine_speed .* cos(fine_phi), fine_speed .* sin(fine_phi), 'b-', ...
    'LineWidth', 1.3, 'DisplayName', 'Interpolated');
hold on;
plot(speed .* cos(phi), speed .* sin(phi), 'ko', ...
    'MarkerFaceColor', 'k', 'MarkerSize', 5, 'DisplayName', 'Computed points');
grid on; box on;
xlabel('$v_{\infty,x}$ [km/s]', 'Interpreter', 'latex');
ylabel('$v_{\infty,z}$ [km/s]', 'Interpreter', 'latex');
legend('Location', 'west');
xlim([-0.4, 3.4]);
ylim([0, 2.5]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
end

function handle = local_cloud(result, visible)
%LOCAL_CLOUD Render the symmetry-reconstructed speed surface from the polar sweep.

data = result.data;
handle = rpo.create_figure(visible,[900,700]);
valid = data.converged & isfinite(data.v_inf_magnitude_km_s);
if nnz(valid) < 2
    text(0.5, 0.5, 'Insufficient converged directions for symmetry reconstruction', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
    return
end
phi = data.polar_angle_rad(valid);
speed = data.v_inf_magnitude_km_s(valid);
[phi, order] = sort(phi);
speed = speed(order);
[phi, unique_indices] = unique(phi, 'stable');
speed = speed(unique_indices);
n_lambda = result.data.renderer.azimuth_samples;
n_beta = result.data.renderer.polar_samples;
[lambda_grid, beta_grid] = meshgrid(linspace(0, 2 * pi, n_lambda), ...
    linspace(-pi / 2, pi / 2, n_beta));
x = cos(beta_grid) .* cos(lambda_grid);
y = cos(beta_grid) .* sin(lambda_grid);
z = sin(beta_grid);
speed_grid = interp1(phi, speed, asin(abs(z)), 'pchip', 'extrap');
surf(speed_grid .* x, speed_grid .* y, speed_grid .* z, speed_grid, ...
    'EdgeColor', 'none');
axis equal; grid on; box on; view(35, 22);
xlabel('$v_{\infty,x}$ [km/s]', 'Interpreter', 'latex');
ylabel('$v_{\infty,y}$ [km/s]', 'Interpreter', 'latex');
zlabel('$v_{\infty,z}$ [km/s]', 'Interpreter', 'latex');
cb = colorbar;
cb.Label.String = '$|v_{\infty}|$ [km/s]';
cb.Label.Interpreter = 'latex';
colormap(parula);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
end
