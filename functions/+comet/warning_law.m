function result = warning_law(settings)
%WARNING_LAW Construct raw and calibrated piecewise-uniform warning-time laws.
%
% PROTOTYPE:
%   result = comet.warning_law(settings)
%
% DESCRIPTION:
%   Reverse the ten survey percentages into chronological bins,
%   split the second year at the first threshold, and rescale four
%   interval groups to the requested exceedance probabilities. Normalize
%   the calibrated masses and integrate each uniform segment exactly.
%
% INPUT:
%   settings  Comet settings with warning.survey_percentages [1 x 10],
%               thresholds_years [1 x 3] and exceedances [1 x 3].
%               Frozen thresholds are 400/365.25, 2 and 3 years, with
%               exceedances 0.90, 0.80 and 0.75. The grouping is fixed
%               for these thresholds and ten one-year survey bins.
%
% OUTPUT:
%   result    Struct with settings, edges_years [1 x 12], raw_masses
%               and baseline_masses [1 x 11, probability per segment],
%               raw/baseline_yearly_survey [1 x 10, reverse chronological
%               survey order], threshold_table (years and exceedances),
%               yearly_table (chronological bin starts and probabilities)
%               and units metadata. Time is years, probabilities are [-].

cfg = settings.warning;
raw_year = fliplr(cfg.survey_percentages/sum(cfg.survey_percentages));
edges = [0,1,cfg.thresholds_years(1),2:10];
raw = [raw_year(1),raw_year(2)*diff(edges(2:3)), ...
    raw_year(2)*diff(edges(3:4)),raw_year(3:end)];

% Match threshold exceedances while retaining relative masses within each group.
groups = {1:2,3,4,5:11};
target = [1-cfg.exceedances(1),-diff(cfg.exceedances),cfg.exceedances(end)];
baseline = raw;
for k = 1:4
    indices = groups{k};
    baseline(indices) = raw(indices)*target(k)/sum(raw(indices));
end
baseline = baseline/sum(baseline);
baseline_year = [baseline(1),sum(baseline(2:3)),baseline(4:end)];
result = struct('settings',cfg,'edges_years',edges,'raw_masses',raw, ...
    'baseline_masses',baseline,'raw_yearly_survey',fliplr(raw_year), ...
    'baseline_yearly_survey',fliplr(baseline_year));

% Integrate the segment laws at thresholds and restore chronological yearly bins.
fractions = max(0,min(1,(edges(2:end)-cfg.thresholds_years(:))./diff(edges)));
result.threshold_table = table(cfg.thresholds_years(:),fractions*raw(:),fractions*baseline(:), ...
    'VariableNames',{'threshold_years','raw_exceedance','baseline_exceedance'});
result.yearly_table = table((0:9).',raw_year(:),baseline_year(:), ...
    'VariableNames',{'warning_bin_start_years','raw_probability','baseline_probability'});
result.units = struct('time','years','probability','dimensionless');
end
