function figures = plot_access(result, visibility)
%PLOT_ACCESS Recreate the historical TOF-C3 maps and example connections.
%
% PROTOTYPE:
%   figures = earth2rpo.plot_access(result, visibility)
%
% DESCRIPTION:
%   Plot retained 1000/10000 km boundary clouds, redraw low-impulse cases on
%   top, and show two stored connection geometries. Preserve the historical
%   C3 quantity without multiplying by the squared velocity unit.
%
% INPUT:
%   result      Finalized study result supplying data.maps, examples, rpo_reference, and configuration.earth_to_rpo.shared.
%   visibility  Figure visibility, on or off.
%
% OUTPUT:
%   figures     Four registered figure handles: two access maps and two example geometries.
%                 Handles are ready for the shared numbered exporter.

s = result.configuration.earth_to_rpo.shared;
data = result.data;

% Render the two retained boundary clouds using the historical C3/impulse scales.
for altitude = [1000,10000]
    map = data.maps.(sprintf('h%d', altitude));
    name = sprintf('TOFvsC3_%d', altitude);
    if altitude == 1000, name = [name '_highlights']; end
    figures.(name) = rpo.create_figure(visibility,[900,650],'Name',name);
    scatter(-map.time_tu*s.TU_days, map.c3_historical, 50, ...
        map.delta_v_tu/s.TU_days/3600/24*s.DU*1000, 'filled');
    hold on; grid on; box on;
    % The historical plot redraws the small-impulse cases on top of the cloud.
    impulse_mps = map.delta_v_tu/s.TU_days/3600/24*s.DU*1000;
    small_impulse = impulse_mps < 12;
    scatter(-map.time_tu(small_impulse)*s.TU_days,map.c3_historical(small_impulse), ...
        50,impulse_mps(small_impulse),'filled');
    xlabel('Time of Flight (days)', 'Interpreter', 'latex');
    % The archived C3 values have not been multiplied by the velocity unit squared.
    ylabel('C$_3$ (nondimensional)', 'Interpreter', 'latex');
    title(sprintf('TOF vs C$_3$ with Color Coded $\\Delta V$ (10$^%d$ km altitude)', ...
        round(log10(altitude))),'Interpreter','latex');
    colormap(gray(256)); clim([0 60]);
    c = colorbar; c.Ticks = 0:10:60;
    c.Label.String = '$\Delta V$ (m/s)'; c.Label.Interpreter = 'latex';
    xlim([0,220]);
    set(gca,'FontName','Times New Roman','FontSize',18);
    if altitude == 1000
        for k = 1:2
            e = data.examples(k);
            plot(e.tof_days, e.c3_historical, 'rp', 'MarkerSize',15,'LineWidth',1.5);
        end
        text(100,-1.3,'(i)','Color','r','FontSize',20,'FontWeight','bold','Interpreter','latex');
        text(101,-1.8,'(ii)','Color','r','FontSize',20,'Interpreter','latex');
    end
end

% Show the frozen manuscript connections at the full and close-up map scales.
for k = 1:2
    name = sprintf('Sample_solution_%d',k);
    e = data.examples(k);
    figures.(name) = rpo.create_figure(visibility,[1200,550],'Name',name);
    for panel = 1:2
        subplot(1,2,panel); hold on; grid on; box on;
        style = '-'; if panel == 2, style = '--'; end
        plot(data.rpo_reference(:,1),data.rpo_reference(:,2),style,'Color',[.5 .5 .5],'LineWidth',2);
        plot(e.state_rot(:,1),e.state_rot(:,2),'b','LineWidth',1);
        plot([-s.mu,1-s.mu],[0,0],'k.','MarkerSize',15);
        xlabel('X (EMB centered rotating frame)','Interpreter','latex');
        ylabel('Y (EMB centered rotating frame)','Interpreter','latex');
        axis equal;
        if panel == 1
            xlim([-5.5,5.5]); ylim([-5.5,5.5]);
            text(-.5,.7,'RPO','Color',[.5 .5 .5],'Interpreter','latex');
            text(-.5,-.3,'Earth','Interpreter','latex'); text(1.1,0,'Moon','Interpreter','latex');
        else
            xlim([-1.5,1.2]); ylim([-1.2,1.2]);
            text(-.2,.95,'RPO','Color',[.5 .5 .5],'Interpreter','latex');
            text(-.1,-.1,'Earth','Interpreter','latex'); text(1,0,'Moon','Interpreter','latex');
            tip = data.rpo_reference(1,1:2); tail = tip - .3*[cosd(135),sind(135)];
            quiver(tail(1),tail(2),tip(1)-tail(1),tip(2)-tail(2),0,'k','LineWidth',1.5,'MaxHeadSize',1);
            text(tail(1),tail(2),'insertion','FontSize',14,'Interpreter','latex');
        end
        set(gca,'FontName','Times New Roman','FontSize',16);
    end
    labels = {'(i)','(ii)'};
    sgtitle(sprintf('%s TOF = %.3f (days), C_3 = %.5f (nondimensional), dv = %.3f (m/s)', ...
        labels{k},e.tof_days,e.c3_historical,e.impulse_mps),'FontName','Times New Roman','FontSize',18);
end
end
