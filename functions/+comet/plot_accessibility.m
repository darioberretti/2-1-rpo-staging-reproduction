function figures = plot_accessibility(result,visibility)
%PLOT_ACCESSIBILITY Draw four full-range and four limited-time accessibility maps.
%
% PROTOTYPE:
%   figures = comet.plot_accessibility(result, visibility)
%
% DESCRIPTION:
%   Render the supplied result using the fixed Run 06 figure
%   names; rpo2ds.export_figures handles numbering and file export.
%
% INPUT:
%   result      Output of comet.compute_accessibility or saved replay;
%                see that producer for the consumed field schema.
%   visibility  Figure visibility: on or off.
%
% OUTPUT:
%   figures     Struct with registered figure names:
%                earth/l1/l2/l1_l2_full_tof_accessibility and the corresponding
%                *_limited_time_accessibility handles; thresholds use km/s
%                and rotating heliocentric positions use AU.

figures = struct();
stems = ["earth","l1","l2","l1_l2"];
for index = [4,2,3,1]
    map = result.maps(index);
    stem = stems(index)+"_full_tof_accessibility";
    figures.(stem) = rpo.create_figure(visibility,[650,600],'Name',map.label+" accessibility");
    panel(map,map.cost_kms(:,:,1),result.settings.accessibility_thresholds_kms);
    title(map.label+" - Unlimited ToF",'Interpreter','none');
    color_scale(gca,result.settings.accessibility_thresholds_kms,false);
    stem = stems(index)+"_limited_time_accessibility";
    figures.(stem) = rpo.create_figure(visibility,[1200,400],'Name',map.label+" limited-time accessibility");
    layout = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    for k = 1:3
        nexttile(layout);
        panel(map,map.cost_kms(:,:,k+1),result.settings.accessibility_thresholds_kms);
        title(sprintf('$T_{\\max}=%.0f$ days',map.horizons_days(k+1)),'Interpreter','latex');
    end
    color_scale(gca,result.settings.accessibility_thresholds_kms,true);
end
end

function panel(map,cost,thresholds)
%PANEL Draw the two cost bands over rotating heliocentric encounter geometry.
%
% PROTOTYPE:
%   panel(map, cost, thresholds)
%
% DESCRIPTION:
%   Transpose radius/angle cost storage for the Cartesian mesh.
%   Leave costs above the upper threshold or nonfinite values transparent.
%
% INPUT:
%   map         Accessibility map with radius_au [R x 1] and
%                angle_deg [A x 1] grids from compute_accessibility.
%   cost        [R x A] Required excess speed [km/s].
%   thresholds  [1 x 2] Cost-band upper limits [km/s].
%
% OUTPUT:
%   None. Draws the cost bands and Sun/Earth/radius guides on gca [AU].

[radius,angle] = meshgrid(map.radius_au,map.angle_deg);
values = NaN(size(cost));
values(cost<=thresholds(1)) = thresholds(1)/2;
values(cost>thresholds(1)&cost<=thresholds(2)) = mean(thresholds);
handle = surf(radius.*cosd(angle),radius.*sind(angle),zeros(size(radius)),values.', ...
    'EdgeColor','none','FaceColor','flat');
set(handle,'AlphaData',double(isfinite(values.')),'FaceAlpha','flat','AlphaDataMapping','none');
view(2); hold on;
theta = linspace(0,2*pi,721);
for radius_value = [map.radius_au(1),1,map.radius_au(end)]
    color = [0.25,0.25,0.25];
    if radius_value==1, color=[0,0.35,0.25]; end
    plot(radius_value*cos(theta),radius_value*sin(theta),'Color',color,'LineWidth',0.8);
end
plot([0,map.radius_au(end)],[0,0],'k-','LineWidth',0.8);
plot([0,0],[0,map.radius_au(end)],'k-','LineWidth',0.8);
plot(0,0,'o','MarkerFaceColor',[1,0.78,0.2],'MarkerEdgeColor','none','MarkerSize',18);
plot(1,0,'ko','MarkerFaceColor',[0.2,0.75,1],'MarkerSize',7);
axis equal; grid on; set(gca,'Layer','top','FontSize',12);
xlim([-1.35,1.55]); ylim([-1.35,1.55]);
xlabel('X Coordinate [au]','Interpreter','latex'); ylabel('Y Coordinate [au]','Interpreter','latex');
colormap(gca,[0.62,0.10,0.70;1,0.94,0]); clim([0,thresholds(2)]);
end

function color_scale(axis_handle,thresholds,shared)
%COLOR_SCALE Label cost thresholds [km/s]; shared=true moves the colorbar to the tile edge.

bar_handle = colorbar(axis_handle);
if shared, bar_handle.Layout.Tile='east'; end
bar_handle.Label.String = '$|v_{\infty}^{+}|$ [km/s]';
bar_handle.Label.Interpreter = 'latex';
bar_handle.Ticks = [0,thresholds];
end
