function figures = plot_proxy(result,visibility)
%PLOT_PROXY Draw the sampled encounter cloud and exact magnitude-mismatch frontier.
%
% PROTOTYPE:
%   figures = sel2rpo.plot_proxy(result, visibility)
%
% DESCRIPTION:
%   Render the saved result using the fixed Run 04 figure names.
%   Export and figure numbering belong to rpo2ds.export_figures.
%
% INPUT:
%   result      Result from sel2rpo.compute_proxy or saved replay:
%                data.plot_sample and pareto_front tables with tof_days
%             [days] and vinf_diff_m_s [m/s].
%   visibility  Figure visibility: on or off.
%
% OUTPUT:
%   figures     Struct of figure handles keyed by registered names:
%                sel2_rpo_database_proxy_pareto; this is a scalar proxy plot.

data = result.data;
fig = rpo.create_figure(visibility,[850,600]);
ax = axes(fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
scatter(ax,data.plot_sample.tof_days,data.plot_sample.vinf_diff_m_s, ...
    7,[0.7,0.7,0.7],'filled','MarkerFaceAlpha',0.18,'MarkerEdgeAlpha',0.18);
plot(ax,data.pareto_front.tof_days,data.pareto_front.vinf_diff_m_s,'r-o', ...
    'LineWidth',1.4,'MarkerSize',4.5,'MarkerFaceColor',[0.85,0.10,0.10],'MarkerEdgeColor','none');
xlabel(ax,'Transfer TOF lower bound [days]');
ylabel(ax,'$\left|\|\mathbf{v}_{\infty,M}^{-}\|-\|\mathbf{v}_{\infty,M}^{+}\|\right|$ [m/s]', ...
    'Interpreter','latex');
set(ax,'FontName','Times New Roman','FontSize',12);
legend(ax,{'Database-pair sample','Proxy Pareto front'},'Location','northeast');
figures.sel2_rpo_database_proxy_pareto = fig;
end
