function figures = plot_validation(result, visibility)
%PLOT_VALIDATION Plot the seed, shooting nodes, and optimized transfer views.
%
% PROTOTYPE:
%   figures = earth2rpo.plot_validation(result, visibility)
%
% DESCRIPTION:
%   Create patch overviews and stored node/solution views at the frozen zoom
%   levels. Rotating views use DU; the inertial plotting helper also uses
%   normalized positions. Display-only arrow scales do not alter stored burns.
%
% INPUT:
%   result      Finalized study result supplying data.seed, initial_setup/initial_guess, setup/solution,
%                  and configuration.earth_to_rpo.patch.
%   visibility  Figure visibility, on or off.
%
% OUTPUT:
%   figures     Eleven registered figure handles for patch, node, rotating/inertial views.
%                 Handles are ready for the shared numbered exporter.

data = result.data;
seed = data.seed;

% Show the selected two-body/database patch at full and Earth-centered scales.
for zoom = [false,true]
    name = 'initial_guess_full';
    if zoom, name = 'initial_guess_zoom'; end
    figures.(name) = rpo.create_figure(visibility,[900,600],'Name',name);
    draw_patch(seed.plot_data);
    axis equal;
    if zoom
        s = seed.constants;
        width = result.configuration.earth_to_rpo.patch.zoom_half_width_km/s.DU;
        xlim(-s.mu + [-width(1),width(1)]); ylim([-width(2),width(2)]);
    else
        p = seed.plot_data;
        x = [p.xx_rpo(:,1);p.vleo_x(:);p.xx_db(:,1);p.xx_conic(:,1);p.moon_xy(1)];
        y = [p.xx_rpo(:,2);p.vleo_y(:);p.xx_db(:,2);p.xx_conic(:,2);p.moon_xy(2)];
        margin = result.configuration.earth_to_rpo.patch.full_plot_margin_du;
        xlim([min(x)-margin,max(x)+margin]); ylim([min(y)-margin,max(y)+margin]);
    end
end

% Draw the initial shooting layout using the shared stored-solution renderer.
node_names = {'initial_guess_full_nodes','initial_guess_zoom_nodes','initial_guess_zoomzoom_nodes'};
node_limits = {[],[-1.3,2,-1.2,1.2],[-.04,.1,-.04,.04]};
for k = 1:3
    f = rpo.create_figure(visibility,[900,600],'Name',node_names{k});
    plot_vleo_rpo_optimization_2d(data.initial_setup,data.initial_guess, ...
        f.Number,'Initial guess',true,node_limits{k},'',false);
    legend('Location','northeast');
    figures.(node_names{k}) = f;
end

% Show optimized rotating and inertial geometries at the three fixed zoom levels.
for inertial = [false,true]
    frame = 'rot'; if inertial, frame = 'inertial'; end
    limits = {[],[-.04,.1,-.04,.04],[-1.3,2,-1.2,1.2]};
    if inertial, limits = {[],[-.04,.1,-.04,.04],[-2.5,2,-1.1,2]}; end
    suffixes = {'','_leo_closeup','_rpo_closeup'};
    for k = 1:3
        name = ['vleo_rpo_optimized_' frame suffixes{k}];
        f = rpo.create_figure(visibility,[900,600],'Name',name);
        if inertial
            plot_vleo_rpo_optimization_inertial_2d(data.setup,data.solution, ...
                f.Number,'Optimized solution',limits{k},'',false);
        else
            plot_vleo_rpo_optimization_2d(data.setup,data.solution, ...
                f.Number,'Optimized solution',false,limits{k},'',false);
        end
        location = 'northeast';
        if k == 2, location = 'southeast'; end
        legend('Location',location);
        figures.(name) = f;
    end
end
end

function draw_patch(p)
%DRAW_PATCH Draw the stored rotating-frame patch geometry.
%
% DESCRIPTION:
%   p is seed.plot_data; positions use DU, with Earth at x = -p.mu.

hold on; grid on; box on;
plot(p.xx_rpo(:,1),p.xx_rpo(:,2),'--','Color',[.25 .45 .90], ...
    'LineWidth',1.4,'DisplayName','RPO reference');
plot(p.vleo_x,p.vleo_y,'Color',[.05 .55 .45],'LineWidth',1.2,'DisplayName','Departure orbit');
plot(p.xx_db(:,1),p.xx_db(:,2),'Color',[.08 .08 .10],'LineWidth',1.5,'DisplayName','RPO-side arc');
plot(p.xx_conic(:,1),p.xx_conic(:,2),'Color',[.88 .30 .12],'LineWidth',1.5,'DisplayName','Earth conic arc');
plot(-p.mu,0,'.','Color',[.05 .25 .85],'MarkerSize',26,'DisplayName','Earth');
plot(p.moon_xy(1),p.moon_xy(2),'.','Color',[.5 .5 .5],'MarkerSize',22,'DisplayName','Moon');
xlabel('x [DU]'); ylabel('y [DU]'); legend('Location','best','Box','on');
set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',.9, ...
    'GridAlpha',.18,'MinorGridAlpha',.10,'Layer','top');
end
