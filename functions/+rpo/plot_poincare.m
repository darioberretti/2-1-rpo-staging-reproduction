function figures = plot_poincare(result, visible)
%PLOT_POINCARE Plot the positive-crossing phase-space section.
%
% PROTOTYPE:
%   figures = rpo.plot_poincare(result, visible)
%
% DESCRIPTION:
%   Plot retained (x, vx) crossings and mark the reference RPO state.
%
% INPUT:
%   result   Finalized study result supplying
%              data.crossing_x / crossing_vx / reference_state.
%   visible  Figure visibility, on or off.
%
% OUTPUT:
%   figures  Struct of registered figure names and handles, ready for
%              rpo.export_figures. Plotting does not write files.

% Show retained crossings with the reference RPO marked in the same section.
figures.rpo_poincare_section = rpo.create_figure(visible,[]);
plot(result.data.crossing_x(:), result.data.crossing_vx(:), '.', ...
    'Color', [0.7, 0.7, 0.7], 'MarkerSize', 8);
hold on;
plot(result.data.reference_state(1), result.data.reference_state(3), ...
    'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
xlabel('$x$ [DU]', 'Interpreter', 'latex');
ylabel('$v_x$ [VU]', 'Interpreter', 'latex');
title('Phase-space structure around the 2:1 RPO', 'Interpreter', 'latex');
xlim([-1.25, -1.05]);
ylim([-0.012, 0.012]);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
end
