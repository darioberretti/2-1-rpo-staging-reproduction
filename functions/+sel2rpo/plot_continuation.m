function figures = plot_continuation(result,visibility)
%PLOT_CONTINUATION Draw merged delta-v and flight/elapsed times versus Sun phase.
%
% PROTOTYPE:
%   figures = sel2rpo.plot_continuation(result, visibility)
%
% DESCRIPTION:
%   Render the saved result using the fixed Run 04 figure names.
%   Export and figure numbering belong to rpo2ds.export_figures.
%
% INPUT:
%   result      Result from sel2rpo.compute_continuation or saved replay:
%                data.cases with theta_sun0_rad [rad], delta_v_m_s [m/s],
%             flight_days and elapsed_days [days].
%   visibility  Figure visibility: on or off.
%
% OUTPUT:
%   figures     Struct of figure handles keyed by registered names:
%                sel2_rpo_continuation_delta_v and sel2_rpo_continuation_times.

cases = result.data.cases;
theta = [cases.theta_sun0_rad];

% Plot the merged physical impulse cost against wrapped Sun phase.
figures.sel2_rpo_continuation_delta_v = rpo.create_figure(visibility,[900,500]);
ax = axes(figures.sel2_rpo_continuation_delta_v);
plot(ax,theta,[cases.delta_v_m_s],'-ok','LineWidth',1,'MarkerSize',4);
ylabel(ax,'Total $\Delta V$ [m/s]','Interpreter','latex'); style(ax);

% Distinguish flight after departure from total time including halo waiting.
figures.sel2_rpo_continuation_times = rpo.create_figure(visibility,[900,500]);
ax = axes(figures.sel2_rpo_continuation_times); hold(ax,'on');
plot(ax,theta,[cases.flight_days],'-ok','LineWidth',1,'MarkerSize',4);
plot(ax,theta,[cases.elapsed_days],'--','Color',[0.25,0.45,0.75],'LineWidth',1.3);
ylabel(ax,'Time [days]'); style(ax);
legend(ax,{'Flight after departure','Total elapsed, including halo waiting'},'Location','best');
end

function style(ax)
%STYLE Style continuation axes over Sun phase [0,2*pi] rad.

grid(ax,'on'); box(ax,'on'); xlim(ax,[0,2*pi]);
xlabel(ax,'$\theta_{S,0}$ [rad]','Interpreter','latex');
set(ax,'FontName','Times New Roman','FontSize',13);
end
