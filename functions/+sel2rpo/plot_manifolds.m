function figures = plot_manifolds(result,visibility)
%PLOT_MANIFOLDS Draw both unstable-velocity branch illustrations in SE projections.
%
% PROTOTYPE:
%   figures = sel2rpo.plot_manifolds(result, visibility)
%
% DESCRIPTION:
%   Render the saved result using the fixed Run 04 figure names.
%   Export and figure numbering belong to rpo2ds.export_figures.
%
% INPUT:
%   result      Result from sel2rpo.compute_manifolds or saved replay:
%                data.arcs.plus/minus (SE state histories), halo_states
%             [M x 6] and constants.mu_SE_nd [-]; positions use DU_SE.
%   visibility  Figure visibility: on or off.
%
% OUTPUT:
%   figures     Struct of figure handles keyed by registered names:
%                sel2_unstable_direction_projections with x-y and x-z panels.

data = result.data;
fig = rpo.create_figure(visibility,[1100,620]);
coordinates = [2,3];
coordinate_labels = {'$y$ [-]','$z$ [-]'};
panel_titles = {'x-y','x-z'};
positions = [0.095,0.105,0.35,0.84;0.56,0.105,0.375,0.84];
vertical_limits = [-0.0255,0.0165;-0.019,0.0197];
for panel = 1:2
    ax = subplot(1,2,panel,'Parent',fig);
    hold(ax,'on');
    for phase = 1:numel(data.arcs.plus)
        positive = data.arcs.plus{phase}.xx_nd;
        negative = data.arcs.minus{phase}.xx_nd;
        plot(ax,positive(:,1),positive(:,coordinates(panel)),'Color',[0.85,0.2,0.2],'LineWidth',0.5);
        plot(ax,negative(:,1),negative(:,coordinates(panel)),'Color',[0.2,0.35,0.85],'LineWidth',0.5);
    end
    plot(ax,data.halo_states(:,1),data.halo_states(:,coordinates(panel)),'k-','LineWidth',1.5);
    plot(ax,1-data.constants.mu_SE_nd,0,'o','MarkerFaceColor',[0.2,0.7,0.2], ...
        'MarkerEdgeColor','k','MarkerSize',6);
    grid(ax,'on'); box(ax,'on'); axis(ax,'equal');
    set(ax,'FontName','Times New Roman','FontSize',12);
    xlabel(ax,'$x$ [-]','Interpreter','latex');
    ylabel(ax,coordinate_labels{panel},'Interpreter','latex');
    title(ax,panel_titles{panel},'FontWeight','normal','FontSize',13);
    set(ax,'Position',positions(panel,:));
    xlim(ax,[0.9925,1.0235]); ylim(ax,vertical_limits(panel,:));
end
figures.sel2_unstable_direction_projections = fig;
end
