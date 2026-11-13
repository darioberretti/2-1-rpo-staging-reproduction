function figures = plot_continuation(result, visibility)
%PLOT_CONTINUATION Plot the stored departure/insertion maneuver breakdown.
%
% PROTOTYPE:
%   figures = earth2rpo.plot_continuation(result, visibility)
%
% DESCRIPTION:
%   Draw a stacked bar per stored altitude and label total maneuver cost.
%   This displays the saved campaign without rerunning a solve.
%
% INPUT:
%   result      Finalized study result supplying data.results_table with altitude_km and departure/insertion/total costs.
%   visibility  Figure visibility, on or off.
%
% OUTPUT:
%   figures     One registered altitude-sweep figure handle.
%                 Handles are ready for the shared numbered exporter.

% Stack departure and insertion maneuvers in the stored altitude order.
t = result.data.results_table;
figures.leo_rpo_altitude_sweep_stacked_dv = rpo.create_figure(visibility, ...
    [900,550],'Name','Departure-altitude continuation');
b = bar(1:height(t),[t.departure_dv_mps,t.injection_dv_mps],'stacked');
b(1).FaceColor = [.20 .45 .85]; b(2).FaceColor = [.90 .45 .15];
grid on; box on;
xticks(1:height(t)); xticklabels(string(t.altitude_km));
xlabel('Departure orbit altitude [km]','Interpreter','latex');
ylabel('$\Delta V$ [m/s]','Interpreter','latex');
legend({'Earth departure maneuver','RPO insertion maneuver'},'Location','northwest');
set(gca,'FontName','Times New Roman','FontSize',13);

% Label total costs with a display offset derived from the plotted range.
offset = .02*max(t.total_dv_mps);
for k = 1:height(t)
    text(k,t.total_dv_mps(k)+offset,sprintf('%.0f',t.total_dv_mps(k)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontName','Times New Roman','FontSize',10);
end
ylim([2500,max(t.total_dv_mps)+4*offset]);
end
