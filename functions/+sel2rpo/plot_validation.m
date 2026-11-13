function figures = plot_validation(result,visibility)
%PLOT_VALIDATION Draw initial and optimized shooting arcs in EM and SE frames.
%
% PROTOTYPE:
%   figures = sel2rpo.plot_validation(result, visibility)
%
% DESCRIPTION:
%   Render the saved result using the fixed Run 04 figure names.
%   Export and figure numbering belong to rpo2ds.export_figures.
%
% INPUT:
%   result      Result from sel2rpo.compute_validation or saved replay:
%                data.initial_guess, solution and problem; trajectories
%             follow sel2rpo.reconstruct_transfer (normalized state rows).
%   visibility  Figure visibility: on or off.
%
% OUTPUT:
%   figures     Struct of figure handles keyed by registered names:
%                sel2_rpo_initial_guess, sel2_rpo_multiple_shooting_nodes,
%                sel2_rpo_optimized_sun_earth and
%                sel2_rpo_optimized_earth_moon.

% Keep initial-guess discontinuities visible in both construction views.
data = result.data;
figures.sel2_rpo_initial_guess = transfer_view(data.initial_guess,data.problem,'em',false,visibility);
title(gca,'Reconstructed initial guess','FontWeight','normal');
figures.sel2_rpo_multiple_shooting_nodes = transfer_view(data.initial_guess,data.problem,'em',true,visibility);
title(gca,'Reconstructed initial guess: shooting nodes','FontWeight','normal');

% Show the optimized arcs and nodes in both normalized rotating frames.
figures.sel2_rpo_optimized_sun_earth = transfer_view(data.solution,data.problem,'se',true,visibility);
figures.sel2_rpo_optimized_earth_moon = transfer_view(data.solution,data.problem,'em',true,visibility);
end

function fig = transfer_view(trajectory,problem,frame,show_nodes,visibility)
%TRANSFER_VIEW Draw separate shooting arcs and optional nodes in a selected rotating frame.
%
% PROTOTYPE:
%   fig = transfer_view(trajectory, problem, frame, show_nodes, visibility)
%
% DESCRIPTION:
%   Plot each arc without connecting discontinuous endpoints.
%   In the SE view, add the circular lunar-orbit guide about the
%   historical secondary center; EM views mark Earth and Moon.
%
% INPUT:
%   trajectory  Reconstructed transfer; see sel2rpo.reconstruct_transfer.
%   problem     Frozen problem with nlp.mu_em, mu_se_nd and length scales.
%   frame       em or se; plotted x-y positions use DU_EM or DU_SE.
%   show_nodes  Logical controlling shooting-node markers.
%   visibility  Figure visibility: on or off.
%
% OUTPUT:
%   fig         Figure handle containing the selected x-y projection.

fig = rpo.create_figure(visibility,[780,700]);
ax = axes(fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
for segment = 1:numel(trajectory.segments)
    if strcmp(frame,'se')
        states = trajectory.segments{segment}.states_se;
    else
        states = trajectory.segments{segment}.states_em;
    end
    plot(ax,states(:,1),states(:,2),'k-','LineWidth',1);
end
if strcmp(frame,'se')
    nodes = trajectory.nodes_se;
    angle = linspace(0,2*pi,500);
    radius = (1-problem.nlp.mu_em)*problem.nlp.l_em_km/problem.nlp.l_se_km;
    plot(ax,1-problem.nlp.mu_se_nd+radius*cos(angle),radius*sin(angle), ...
        '--','Color',[0.5,0.5,0.5],'LineWidth',0.8);
    plot(ax,1-problem.nlp.mu_se_nd,0,'.b','MarkerSize',12);
    xlabel(ax,'$x$ [$DU_{SE}$]','Interpreter','latex');
    ylabel(ax,'$y$ [$DU_{SE}$]','Interpreter','latex');
    axis(ax,'equal'); xlim(ax,[0.9901,1.0111]); ylim(ax,[-0.0115,0.0075]);
else
    nodes = trajectory.nodes_em;
    plot(ax,-problem.nlp.mu_em,0,'.b','MarkerSize',14);
    plot(ax,1-problem.nlp.mu_em,0,'.','Color',[0.6,0.6,0.6],'MarkerSize',12);
    xlabel(ax,'$x$ [$DU_{EM}$]','Interpreter','latex');
    ylabel(ax,'$y$ [$DU_{EM}$]','Interpreter','latex');
    axis(ax,'equal');
end
if show_nodes
    plot(ax,nodes(:,1),nodes(:,2),'ok','MarkerSize',4.5,'LineWidth',0.8);
end
set(ax,'FontName','Times New Roman','FontSize',13);
end
