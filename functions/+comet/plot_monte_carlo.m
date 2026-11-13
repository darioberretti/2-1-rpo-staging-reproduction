function figures = plot_monte_carlo(result,visibility)
%PLOT_MONTE_CARLO Draw completion curves and mission-capped intercept-epoch bars.
%
% PROTOTYPE:
%   figures = comet.plot_monte_carlo(result, visibility)
%
% DESCRIPTION:
%   Render the supplied result using the fixed Run 06 figure
%   names; rpo2ds.export_figures handles numbering and file export.
%
% INPUT:
%   result      Output of comet.compute_monte_carlo or saved replay;
%                see that producer for the consumed field schema.
%   visibility  Figure visibility: on or off.
%
% OUTPUT:
%   figures     Struct with registered figure names:
%                monte_carlo_l1_l2 and monte_carlo_mixed figure handles.
%                Intercept epochs use years; probabilities include failures.

figures.monte_carlo_l1_l2 = comparison(result.combined,result.settings,visibility,'L1/L2 budget comparison');
figures.monte_carlo_mixed = comparison(result.mixed,result.settings,visibility,'Paper architecture comparison');
end

function handle = comparison(study,settings,visibility,name)
%COMPARISON Draw an empirical completion curve and capped epoch summaries.
%
% PROTOTYPE:
%   handle = comparison(study, settings, visibility, name)
%
% DESCRIPTION:
%   At integer-year support, divide cumulative successes by all
%   scenarios. Plot stored median/P95 values, including horizon caps.
%
% INPUT:
%   study       combined or mixed study from compute_monte_carlo;
%                durations_years [S x B] retains failed scenarios as NaN,
%                labels identifies columns, table stores epoch quantiles.
%   settings    Monte Carlo settings: n_scenarios=S and horizon [years].
%   visibility  Figure visibility: on/off.
%   name        Figure display name.
%
% OUTPUT:
%   handle      Figure with completion probability and epoch [years] panels.

handle = rpo.create_figure(visibility,[1100,650],'Name',name);
layout = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile(layout); hold on;
styles = {'-','--','-.',':'};
markers = {'o','s','d','^'};
colors = lines(width(study.durations_years));
labels = string(study.labels);

% Keep failed scenarios in the empirical completion-probability denominator.
times = 0:settings.mission_duration_years;
for k = 1:width(study.durations_years)
    probabilities = arrayfun(@(t) sum(study.durations_years(:,k)<=t)/settings.n_scenarios,times);
    plot(times,probabilities,'LineStyle',styles{k},'Marker',markers{k},'Color',colors(k,:), ...
        'LineWidth',1.8,'MarkerSize',8,'DisplayName',labels(k));
end
grid on; box on; xlim([0,settings.mission_duration_years]); ylim([0,1]);
xlabel('Intercept Epoch from RPO injection [years]'); ylabel('Completion Probability');
title('Completion probability by intercept time');
legend('Location','southeast','Interpreter','latex','FontSize',10);
nexttile(layout);

% Display stored unconditional epoch quantiles, capped at the study horizon.
bars = bar([study.table.median_years,study.table.p95_years],'grouped');
bars(1).FaceColor = [0.45,0.45,0.45];
bars(2).FaceColor = [0.86,0.93,1]; bars(2).EdgeColor = [0.35,0.45,0.55];
grid on; box on; ylim([0,settings.mission_duration_years]);
xticklabels(labels); xtickangle(30); set(gca,'TickLabelInterpreter','latex');
ylabel('Intercept Epoch from RPO injection [years]'); title('Mission duration statistics');
legend('Median','95 percent','Location','northwest');
set(findall(handle,'Type','axes'),'FontSize',12);
end
