function result = compute_monte_carlo(settings,directory)
%COMPUTE_MONTE_CARLO Evaluate both frozen architecture comparisons using shared random draws.
%
% PROTOTYPE:
%   result = comet.compute_monte_carlo(settings, directory)
%
% DESCRIPTION:
%   Seed one opportunity catalogue and calibrated warning sample.
%   Available flight time is min(warning, opportunity epoch). Interpolate
%   cumulative required excess speed with no extrapolation for the
%   combined L1/L2, L2 and Earth databases. The mixed study uses L2
%   for its baseline and Earth data for both later architecture limits.
%   The published scalar RPO limit is independent of recalibration.
%
% INPUT:
%   settings  Comet-study cfg.comet_case_study struct returned by
%               reproduction_config(); see config/reproduction_config.m.
%               Required substructures are identified below.
%   settings.monte_carlo defines seed, scenario count, horizon [years],
%       opportunity rate [1/year], radius [AU]/angle [deg] ranges and
%       combined/mixed excess-speed limits [km/s]. constants.year_days
%       converts sampled years to available flight days.
%   directory  Folder containing the heliocentric database products.
%
% OUTPUT:
%   result    Struct with settings, catalogue (comet_build_base_catalog),
%               sampled_warning_years/available_tof_days [E x 1],
%               required_vinf_kms [E x 3, km/s] in combined/L2/Earth order,
%               database identities, setup_table, units and elapsed_seconds.
%               combined/mixed store durations_years and selected_encounters
%               [S x B], table (success fraction, median/P95 epochs and
%               p95_reached), and labels. E is opportunities, S scenarios
%               and B budget columns. Failures remain NaN in duration/
%               encounter arrays and remain in the probability denominator.
%               Unreached median/P95 epochs are capped at the mission horizon.

timer = tic;
cfg = settings.monte_carlo;
warning = comet.warning_law(settings);
law = struct('segment_edges_years',warning.edges_years,'segment_masses',warning.baseline_masses);

% Draw the shared catalogue first, then sample warnings in historical order.
rng(cfg.rng_seed);
catalogue = comet_build_base_catalog(cfg);
sampled_warning_years = comet_sample_warning_times(law,catalogue.total_encounters);
available_days = min(sampled_warning_years,catalogue.opportunity_epoch_years)*settings.constants.year_days;
roles = [4,3,1]; % Combined, L2, Earth; every database sees identical draws.
required = NaN(catalogue.total_encounters,3);
identities = cell(1,3);

% Interpolate each database at the same encounter geometry and available time.
for k = 1:3
    database = comet.read_database(directory,roles(k));
    interpolant = griddedInterpolant({database.Rc_grid_au,database.theta_grid_deg,database.tof_grid_days}, ...
        database.cumulative_best_vinf_mag_kms,'linear','none');
    required(:,k) = interpolant(catalogue.Rc_au,catalogue.theta_deg,available_days);
    identities{k} = database.identity;
    fprintf('Monte Carlo database %d/3 evaluated (%d opportunities).\n',k,catalogue.total_encounters);
    clear database interpolant
end

% Apply the two fixed database/budget mappings without using the calibration sweep.
combined = summarize(catalogue,required(:,ones(1,4)),cfg.combined_budgets_kms,cfg);
mixed = summarize(catalogue,required(:,[2,3,3]),cfg.mixed_budgets_kms,cfg);
combined.labels = compose('$\\Delta v=%.2f$ km/s',cfg.combined_budgets_kms);
mixed.labels = ["L2 baseline ($0.75$ km/s)","L2 + Moon ($1.50$ km/s)","RPO ($4.08$ km/s)"];
result = struct('settings',cfg,'catalogue',catalogue, ...
    'sampled_warning_years',sampled_warning_years,'available_tof_days',available_days, ...
    'required_vinf_kms',required,'databases',[identities{:}],'combined',combined, ...
    'mixed',mixed,'elapsed_seconds',toc(timer));
result.units = struct('time','years','available_tof','days','radius','AU', ...
    'angle','degrees','excess_speed','km/s');
result.setup_table = table(cfg.rng_seed,cfg.n_scenarios,cfg.mission_duration_years, ...
    cfg.mean_rate_per_year,catalogue.total_encounters, ...
    'VariableNames',{'rng_seed','scenarios','mission_years','opportunity_rate_per_year','sampled_opportunities'});
disp(combined.table); disp(mixed.table);
end

function study = summarize(catalogue,requirements,budgets,cfg)
%SUMMARIZE Select each scenario's earliest feasible opportunity for each budget.
%
% PROTOTYPE:
%   study = summarize(catalogue, requirements, budgets, cfg)
%
% DESCRIPTION:
%   Order feasible rows by scenario, intercept epoch and required
%   speed. Use all scenarios in success probabilities and order-statistic
%   median/P95 epochs; unreached quantiles use the mission horizon.
%   These are absolute intercept epochs, rather than flight durations.
%
% INPUT:
%   catalogue     Opportunity catalogue from comet_build_base_catalog.
%   requirements  [E x B] Required excess speeds [km/s], NaN if unavailable.
%   budgets       [1 x B] Accepted excess-speed limits [km/s].
%   cfg           Monte Carlo settings with n_scenarios and
%                   mission_duration_years [years].
%
% OUTPUT:
%   study         durations_years [S x B, years] and selected_encounters
%                   [S x B, catalogue row indices], NaN on failure;
%                   table with limits, success_fraction, median_years,
%                   p95_years and p95_reached. E is opportunities and
%                   S is n_scenarios; failed scenarios remain in S.

durations = NaN(cfg.n_scenarios,numel(budgets));
selected_encounters = NaN(size(durations));
success_fraction = zeros(numel(budgets),1);
median_years = zeros(numel(budgets),1);
p95_years = zeros(numel(budgets),1);
p95_reached = false(numel(budgets),1);

% Select earliest successful encounters and retain failed scenarios as NaN.
for k = 1:numel(budgets)
    valid = isfinite(requirements(:,k)) & requirements(:,k)<=budgets(k);
    candidates = sortrows([catalogue.scenario_id(valid),catalogue.opportunity_epoch_years(valid), ...
        requirements(valid,k),find(valid)],[1,2,3]);
    first = [true;diff(candidates(:,1))~=0];
    selected = candidates(first,:);
    times = NaN(cfg.n_scenarios,1);
    times(selected(:,1)) = selected(:,2);
    selected_encounters(selected(:,1),k) = selected(:,4);
    durations(:,k) = times;
    completed = sort(times(isfinite(times)));
    success_fraction(k) = numel(completed)/cfg.n_scenarios;
    median_years(k) = cfg.mission_duration_years;
    p95_years(k) = cfg.mission_duration_years;
    if numel(completed)>=ceil(0.5*cfg.n_scenarios)
        median_years(k) = completed(ceil(0.5*cfg.n_scenarios));
    end
    if numel(completed)>=ceil(0.95*cfg.n_scenarios)
        p95_years(k) = completed(ceil(0.95*cfg.n_scenarios));
        p95_reached(k) = true;
    end
end
study = struct('durations_years',durations,'selected_encounters',selected_encounters,'table', ...
    table(budgets(:),success_fraction,median_years,p95_years,p95_reached, ...
    'VariableNames',{'excess_speed_limit_kms','success_fraction','median_years','p95_years','p95_reached'}));
end
