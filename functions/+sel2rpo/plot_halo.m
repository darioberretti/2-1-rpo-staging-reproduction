function figures = plot_halo(result,visibility)
%PLOT_HALO Draw the phase-colored SEL2 halo and three projected shadows.
%
% PROTOTYPE:
%   figures = sel2rpo.plot_halo(result, visibility)
%
% DESCRIPTION:
%   Render the saved result using the fixed Run 04 figure names.
%   Export and figure numbering belong to rpo2ds.export_figures.
%
% INPUT:
%   result      Result from sel2rpo.compute_halo or saved replay:
%                data.trajectory.xx_nd [M x 6, SE rotating DU_SE/VU_SE],
%             t_days [M x 1, days] and constants.mu_SE_nd [-].
%   visibility  Figure visibility: on or off.
%
% OUTPUT:
%   figures     Struct of figure handles keyed by registered names:
%                sel2_reference_halo; color encodes days along the halo.

halo = result.data;
states = halo.trajectory.xx_nd;
phase_days = halo.trajectory.t_days;

% Locate the secondary/L2 markers and pad the three projection planes.
earth = [1-halo.constants.mu_SE_nd,0,0];
[l2_x,~] = libration_point_collinear_cr3bp(halo.constants.mu_SE_nd,'L2');
l2 = [l2_x,0,0];
points = [states(:,1:3);earth;l2];
lower = min(points); upper = max(points);
span = upper-lower;
lower = lower-0.12*span; upper = upper+0.12*span;

% Draw projected shadows and color the spatial halo by elapsed phase days.
fig = rpo.create_figure(visibility,[980,720]);
ax = axes(fig,'Position',[0.075,0.105,0.820,0.815]);
hold(ax,'on');
plot3(ax,states(:,1),states(:,2),repmat(lower(3),size(states,1),1), ...
    'Color',[0.60,0.60,0.58],'LineWidth',1.15);
plot3(ax,states(:,1),repmat(upper(2),size(states,1),1),states(:,3), ...
    'Color',[0.60,0.60,0.58],'LineWidth',1.15);
plot3(ax,repmat(lower(1),size(states,1),1),states(:,2),states(:,3), ...
    'Color',[0.60,0.60,0.58],'LineWidth',1.15);
surface(ax,[states(:,1),states(:,1)],[states(:,2),states(:,2)], ...
    [states(:,3),states(:,3)],[phase_days,phase_days], ...
    'FaceColor','none','EdgeColor','interp','LineWidth',1.45);
colormap(ax,turbo(256)); clim(ax,[0,phase_days(end)]);

% Mark the initial state, secondary and L2 together with their projections.
start = states(1,1:3);
start_handle = scatter3(ax,start(1),start(2),start(3),48,[0,0.48,0.34], ...
    'MarkerFaceColor','w','LineWidth',1.35);
for plane = 1:3
    projected = start;
    if plane==1, projected(3)=lower(3); elseif plane==2, projected(2)=upper(2); else, projected(1)=lower(1); end
    scatter3(ax,projected(1),projected(2),projected(3),34,[0,0.48,0.34], ...
        'MarkerFaceColor','w','LineWidth',1);
end
earth_projection = projections(ax,earth,lower,upper,[0,0.32,0.95],'o',34);
l2_projection = projections(ax,l2,lower,upper,[0.90,0.05,0.04],'s',26);
earth_handle = scatter3(ax,earth(1),0,0,42,[0,0.32,0.95],'filled','MarkerEdgeColor',[0.2,0.2,0.2]);
l2_handle = scatter3(ax,l2(1),0,0,23,[0.90,0.05,0.04],'s','filled');
text(ax,earth(1),0.07*span(2),0.06*span(3),'Earth', ...
    'Interpreter','latex','FontSize',13,'HorizontalAlignment','right');
text(ax,l2(1)+0.02*span(1),-0.06*span(2),0.085*span(3),'$L_2$', ...
    'Interpreter','latex','FontSize',13);
xlim(ax,[lower(1),upper(1)]); ylim(ax,[lower(2),upper(2)]); zlim(ax,[lower(3),upper(3)]);

% Apply the saved viewing geometry, physical-time color scale and legend.
grid(ax,'on'); box(ax,'on'); view(ax,37,18); camproj(ax,'perspective');
pbaspect(ax,[1.45,1,0.95]);
set(ax,'FontName','Times New Roman','FontSize',14,'GridAlpha',0.22);
xlabel(ax,'$x_{SE}$ [-]','Interpreter','latex');
ylabel(ax,'$y_{SE}$ [-]','Interpreter','latex');
zlabel(ax,'$z_{SE}$ [-]','Interpreter','latex');
title(ax,'Reference Sun--Earth $L_2$ Halo Orbit','Interpreter','latex','FontSize',12);
color_scale = colorbar(ax); color_scale.Label.String = 'Halo phase [d]';
color_scale.Label.Interpreter = 'latex'; color_scale.TickLabelInterpreter = 'latex';
legend(ax,[start_handle,earth_handle,l2_handle,earth_projection,l2_projection], ...
    {'Initial state','Earth','$L_2$','Earth projection','$L_2$ projection'}, ...
    'Interpreter','latex','Location','northeast','FontSize',11);
figures.sel2_reference_halo = fig;
end

function handle = projections(ax,point,lower,upper,color,marker,marker_size)
%PROJECTIONS Mark a 3D point on three coordinate planes and return the first marker.
% point, lower and upper are [1 x 3] SE positions [DU_SE]; ax is the axes.

handle = scatter3(ax,point(1),point(2),lower(3),marker_size,color,marker,'LineWidth',1.05);
scatter3(ax,point(1),upper(2),point(3),marker_size,color,marker,'LineWidth',1.05);
scatter3(ax,lower(1),point(2),point(3),marker_size,color,marker,'LineWidth',1.05);
end
