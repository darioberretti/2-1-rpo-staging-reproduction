function base_catalog = comet_build_base_catalog(settings)
%COMET_BUILD_BASE_CATALOG Draw the shared synthetic encounter catalogue in historical RNG order.
%
% PROTOTYPE:
%   base_catalog = comet_build_base_catalog(settings)
%
% DESCRIPTION:
%   Draw Poisson encounter counts first. For each scenario, sample
%   sorted uniform opportunity epochs, then independent uniform radii
%   and angles. Radius itself is uniform over the configured interval,
%   rather than area-uniform in an annulus. The caller seeds the RNG.
%
% INPUT:
%   settings  Monte Carlo settings with mean_rate_per_year [1/year],
%               mission_duration_years [years], n_scenarios, Rc_range_au
%               [1 x 2, AU] and theta_range_deg [1 x 2, degrees].
%
% OUTPUT:
%   base_catalog  Struct with encounter_counts [S x 1], total_encounters
%                   E, and [E x 1] arrays scenario_id, opportunity_epoch_years
%                   [years], Rc_au [AU], theta_deg [deg]. S is n_scenarios.
%                   discovery_epoch_years is a legacy alias of opportunity
%                   epochs; actual sampled warnings are handled separately.

    lambda_total = settings.mean_rate_per_year * settings.mission_duration_years; % expected encounters per scenario
    encounter_counts = sample_poisson_knuth(lambda_total, settings.n_scenarios);

    % Compute the total number of synthetic encounters in the campaign.
    total_encounters = sum(encounter_counts);                           % total number of synthetic encounters

    % Preallocate the catalog arrays.
    scenario_id = zeros(total_encounters, 1);                          % scenario identifier per encounter
    opportunity_epoch_years = zeros(total_encounters, 1);              % intercept opportunity epoch [years]
    Rc_au = zeros(total_encounters, 1);                                % intercept radius [AU]
    theta_deg = zeros(total_encounters, 1);                            % intercept phase angle [deg]

    % Fill the common encounter geometry for all laws.
    row0 = 0;
    for i_scenario = 1:settings.n_scenarios
        n_here = encounter_counts(i_scenario);                         % encounter count in this scenario
        idx = row0 + (1:n_here);                                       % rows of this scenario inside the flat arrays

        % Store the scenario index of the current encounter block.
        scenario_id(idx) = i_scenario;

        % Sample and sort the intercept opportunity epochs inside the mission window.
        opportunity_epoch_years(idx) = sort(settings.mission_duration_years * rand(n_here, 1));

        % Sample the intercept radius uniformly in the allowed annulus.
        Rc_au(idx) = settings.Rc_range_au(1) + diff(settings.Rc_range_au) * rand(n_here, 1);

        % Sample the encounter phase angle uniformly in the allowed interval.
        theta_deg(idx) = settings.theta_range_deg(1) + diff(settings.theta_range_deg) * rand(n_here, 1);

        % Advance the flat-row pointer to the next scenario block.
        row0 = row0 + n_here;
    end

    % Package the synthetic encounter catalog.
    base_catalog = struct();
    base_catalog.encounter_counts = encounter_counts;                  % encounter count per scenario
    base_catalog.total_encounters = total_encounters;                  % total encounter count
    base_catalog.scenario_id = scenario_id;                            % scenario identifier per encounter
    base_catalog.opportunity_epoch_years = opportunity_epoch_years;    % intercept opportunity epoch per encounter [years]
    base_catalog.discovery_epoch_years = opportunity_epoch_years;      % legacy name kept for older callers [years]
    base_catalog.Rc_au = Rc_au;                                        % intercept radius per encounter [AU]
    base_catalog.theta_deg = theta_deg;                                % intercept phase angle per encounter [deg]
end

function counts = sample_poisson_knuth(lambda, n_samples)
%SAMPLE_POISSON_KNUTH Draw n_samples Poisson counts using Knuth's serial uniform-product method.
% lambda is the positive expected count; counts is [n_samples x 1].

    counts = zeros(n_samples, 1);
    threshold = exp(-lambda);                                          % stopping threshold for the direct product method

    % Draw one Poisson random value for each mission scenario.
    for i_sample = 1:n_samples
        k = 0;                                                         % product counter
        p = 1.0;                                                       % cumulative product

        % Multiply uniforms until the Poisson stopping condition is met.
        while p > threshold
            k = k + 1;
            p = p * rand();
        end

        % Convert the loop count into the sampled Poisson value.
        counts(i_sample) = k - 1;
    end
end
