function figures = plot_comparison(result,visibility)
%PLOT_COMPARISON Draw the empirical annual-cost CDFs for the left and right solutions.
%
% PROTOTYPE:
%   figures = stationkeeping.plot_comparison(result, visibility)
%
% DESCRIPTION:
%   Use saved completed-trial statistics directly. Draw a staircase
%   starting at probability zero at the minimum cost; omit empty
%   distributions and keep the fixed historical branch/stability labels.
%
% INPUT:
%   result      stationkeeping.compute result or saved replay, with
%                branches.left/right.statistics.cdf_cost_m_s [M x 1, m/s]
%                and cdf_probability [M x 1, dimensionless]. Each branch
%                can have a different valid-trial count M.
%   visibility  Figure visibility: on or off.
%
% OUTPUT:
%   figures     Struct with stationkeeping_cost_cdf (figure handle).
%                rpo2ds.export_figures assigns its fixed Run 05 number.

fig = rpo.create_figure(visibility,[850,600]);
ax = axes(fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
colors = [0,0.4470,0.7410;0.8500,0.3250,0.0980];
branches = {'left','right'};
labels = {'Left (in-plane unstable)','Right (in-plane stable)'};

% Plot empirical steps without fitting, smoothing or recomputing statistics.
for k = 1:2
    stats = result.branches.(branches{k}).statistics;
    if ~isempty(stats.cdf_cost_m_s)
        stairs(ax,[stats.cdf_cost_m_s(1);stats.cdf_cost_m_s], ...
            [0;stats.cdf_probability],'Color',colors(k,:),'LineWidth',1.6, ...
            'DisplayName',labels{k});
    end
end
set(ax,'FontName','Times New Roman','FontSize',12);
xlabel(ax,'Annual stationkeeping cost [m/s]');
ylabel(ax,'Empirical cumulative probability');
ylim(ax,[0,1]);
legend(ax,'Location','southeast');
figures.stationkeeping_cost_cdf = fig;
end
