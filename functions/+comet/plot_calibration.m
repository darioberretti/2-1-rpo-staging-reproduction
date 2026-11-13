function figures = plot_calibration(result,visibility)
%PLOT_CALIBRATION Historical calibration panels and planar RTH polar boundary.
%
% PROTOTYPE:
%   figures = comet.plot_calibration(result, visibility)
%
% DESCRIPTION:
%   Plot all calibration trials in blue, following the manuscript style.
%   Solver status is retained in result and is not encoded in these plots.
%
% INPUT:
%   result      comet.compute_calibration result or saved replay:
%                 calibration.history has vinf_mag_km_s [km/s], dv_m_s
%                 [m/s] and tof_days [days]; sweep has angle_deg [deg]
%                 and reachable_vinf_mag_km_s [km/s]. settings.first
%                 supplies requested_budget_m_s for the dashed reference.
%   visibility  Figure visibility, 'on' or 'off'.
%
% OUTPUT:
%   figures     Struct with rpo_first_direction_calibration and
%                 rpo_rth_polar_boundary handles. The latter uses unit
%                 Earth-Sun RTH directions, with angles from +R toward +T.

% Show every calibration trial with the established uniform-blue styling.
figures.rpo_first_direction_calibration = rpo.create_figure(visibility,[1000,600], ...
    'Name','RPO first-direction calibration');
history = result.calibration.history;
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;
plot(history.vinf_mag_km_s,history.dv_m_s,'o-','Color',[0,0.4470,0.7410], ...
    'LineWidth',1.4);
hold on;
yline(result.settings.first.requested_budget_m_s,'k--','LineWidth',1);
grid on; box on;
xlabel('$|v_{\infty}|$ [km/s]','Interpreter','latex');
ylabel('Needed $\Delta v$ [m/s]','Interpreter','latex');
title('First-direction budget calibration');
nexttile;
plot(history.vinf_mag_km_s,history.tof_days,'s-','Color',[0,0.4470,0.7410], ...
    'LineWidth',1.4);
grid on; box on;
xlabel('$|v_{\infty}|$ [km/s]','Interpreter','latex'); ylabel('TOF [days]','Interpreter','latex');
title('First-direction transfer time');

% Draw the solved planar Earth-Sun RTH boundary separately from calibration.
figures.rpo_rth_polar_boundary = rpo.create_figure(visibility,[650,600],'Name','RPO RTH boundary');
polarplot(deg2rad(result.sweep.angle_deg),result.sweep.reachable_vinf_mag_km_s,'LineWidth',1.6);
title('Planar RTH reachable $|v_{\infty}|$','Interpreter','latex');
end
