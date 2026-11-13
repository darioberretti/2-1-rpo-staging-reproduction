function figures = plot_warning(result,visibility)
%PLOT_WARNING Draw raw/calibrated warning PMFs, survival laws and thresholds.
%
% PROTOTYPE:
%   figures = comet.plot_warning(result, visibility)
%
% DESCRIPTION:
%   Render the supplied result using the fixed Run 06 figure
%   names; rpo2ds.export_figures handles numbering and file export.
%
% INPUT:
%   result      Output of comet.warning_law or saved replay;
%                see that producer for the consumed field schema.
%   visibility  Figure visibility: on or off.
%
% OUTPUT:
%   figures     Struct with registered figure names:
%                warning_law_generation, a three-panel figure. Time [years]
%                and probabilities [-] follow the result law.

figures.warning_law_generation = rpo.create_figure(visibility,[1300,400], ...
    'Name','Warning-time law');
layout = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
colors = [0.25,0.25,0.25;0,0.447,0.741];
nexttile(layout);
% Each column is one law, with chronological one-year bins on the x-axis.
bars = bar(0:9,[flipud(result.raw_yearly_survey(:)),flipud(result.baseline_yearly_survey(:))]);
for k=1:2, bars(k).FaceColor=colors(k,:); end
grid on; box on;
xlabel('Warning-time bin [yr]'); ylabel('Probability');
title('Yearly warning-time PMF');
legend('Raw Inno law','Baseline calibrated law','Location','northwest');
ylim([0,0.6]);

% Draw exact segment survival functions, followed by threshold exceedances.
nexttile(layout); hold on;
plot(result.edges_years,[1,1-cumsum(result.raw_masses)],'Color',colors(1,:),'LineWidth',1.8);
plot(result.edges_years,[1,1-cumsum(result.baseline_masses)],'Color',colors(2,:),'LineWidth',1.8);
grid on; box on; xlim([0,10]); ylim([0,1.2]);
xlabel('Warning time [yr]'); ylabel('$P(t_{\mathrm{warning}}\geq T)$','Interpreter','latex');
title('Survival functions');
nexttile(layout);
thresholds = result.settings.thresholds_years;
raw = arrayfun(@(t) survival(t,result.edges_years,result.raw_masses),thresholds);
baseline = arrayfun(@(t) survival(t,result.edges_years,result.baseline_masses),thresholds);
bars = bar([raw(:),baseline(:)]);
for k=1:2, bars(k).FaceColor=colors(k,:); end
xticklabels({'400 d','2 yr','3 yr'}); grid on; box on; ylim([0,1.05]);
ylabel('Exceedance probability'); title('Calibration thresholds');
set(findall(figures.warning_law_generation,'Type','axes'),'FontSize',12);
end

function value = survival(t,edges,masses)
%SURVIVAL Integrate a piecewise-uniform law above scalar t [years].
% edges [1 x K+1] are years; masses [1 x K] are segment probabilities.

fraction = max(0,min(1,(edges(2:end)-t)./diff(edges)));
value = sum(masses.*fraction);
end
