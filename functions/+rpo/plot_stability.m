function figures = plot_stability(result, visible)
%PLOT_STABILITY Plot planar and spatial A1 stability indicators.
%
% PROTOTYPE:
%   figures = rpo.plot_stability(result, visible)
%
% DESCRIPTION:
%   Color the period/altitude branches by nontrivial spectral radius,
%   with stable segments and the prescribed resonance marked.
%
% INPUT:
%   result   Finalized study result supplying
%              data branch numbers, periods, altitudes, spectral radii, and configuration.
%   visible  Figure visibility, on or off.
%
% OUTPUT:
%   figures  Struct of registered figure names and handles, ready for
%              rpo.export_figures. Plotting does not write files.

% Use the same period/altitude layout for planar and spatial spectral radii.
figures.a1_planar_stability = local_plot(result, ...
    result.data.planar_spectral_radius, '2D', visible);
figures.a1_spatial_stability = local_plot(result, ...
    result.data.spatial_spectral_radius, '3D', visible);
end

function figure_handle = local_plot(result, spectral_radius, label, visible)
%LOCAL_PLOT Plot one spectral-radius field over the A1 branches.
%
% PROTOTYPE:
%   figure_handle = local_plot(result, spectral_radius, label, visible)
%
% DESCRIPTION:
%   Interpolate each branch in perilune altitude for the color display and
%   overlay its original samples. Stable points use the frozen tolerance;
%   the plot marks the prescribed resonance and the Moon singularity.
%
% INPUT:
%   result           Finalized stability result with data and configuration.
%   spectral_radius  [N x 1] Nontrivial spectral radius for each family row [-].
%   label            2D or 3D title label.
%   visible          Figure visibility, on or off.
%
% OUTPUT:
%   figure_handle    Handle to the period [days] versus altitude [km] figure.

cfg = result.configuration;
settings = cfg.reference_rpo.stability;
figure_handle = rpo.create_figure(visible,[]);
hold on;
grid on;

% Interpolate within each branch and highlight the tolerance-based stable segments.
for branch = 1:3
    mask = result.data.branch_number == branch;
    [altitude, order] = sort(result.data.perilune_altitude_km(mask), 'descend');
    period = result.data.period_days(mask);
    period = period(order);
    radius = spectral_radius(mask);
    radius = radius(order);

    interpolated_altitude = linspace(min(altitude), max(altitude), ...
        settings.interpolation_points);
    interpolated_period = pchip(altitude, period, interpolated_altitude);
    interpolated_radius = pchip(altitude, radius, interpolated_altitude);
    scatter(interpolated_altitude, interpolated_period, 15, ...
        min(interpolated_radius, settings.color_limit), 'filled');
    scatter(altitude, period, 15, min(radius, settings.color_limit), 'filled');

    stable = interpolated_radius <= 1 + settings.tolerance;
    plot(interpolated_altitude(stable), interpolated_period(stable), '.', ...
        'Color', [0.3, 0.9, 1.0], 'MarkerSize', 15);
end

% Mark the prescribed resonance and lunar singularity on the common axes.
colormap(copper);
clim([1, settings.color_limit]);
resonance_days = -pi / cfg.reference_rpo.rpo.synodic_frequency_rad_tu ...
    * cfg.constants.time_unit_days;
yline(resonance_days, 'r--', 'LineWidth', 1.8);
xline(-cfg.constants.moon_radius_km, 'k--', 'LineWidth', 1.2);

% Apply frozen label positions only to the complete 400-orbit reference family.
if numel(result.data.branch_number) == 400
    text(3.55e5, 4, 'Sub-family k1', 'HorizontalAlignment', 'right');
    text(6e4, 48, 'Sub-family k2');
    text(1e5, 57, 'Sub-family k3');
    text(3.55e5, resonance_days + 1.2, '2:1 synodic resonance', ...
        'Color', 'r', 'HorizontalAlignment', 'right');
    text(0, 1.2, 'Moon singularity');
end
colorbar_handle = colorbar;
ylabel(colorbar_handle, '$\max |\lambda|$', 'Interpreter', 'latex');
xlabel('Perilune altitude [km]', 'Interpreter', 'latex');
ylabel('Period [days]', 'Interpreter', 'latex');
title("A1 family linear stability analysis - " + label, ...
    'Interpreter', 'latex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
end
