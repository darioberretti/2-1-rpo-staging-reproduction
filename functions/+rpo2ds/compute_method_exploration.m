function result = compute_method_exploration(cfg, phase_database, altitude_database, targets)
%COMPUTE_METHOD_EXPLORATION Recreate the frozen transfer-construction diagnostics.
%
% PROTOTYPE:
%   result = rpo2ds.compute_method_exploration(cfg, phase_database, ...
%       altitude_database, targets)
%
% DESCRIPTION:
%   Summarize the encounter table, reconstruct illustrative planetary and
%   Earth-departure arcs, and compare two fixed-speed escape directions over
%   perigee altitude. The complete table supplies ranges/counts; only a
%   deterministic subset of its rows is selected for the encounter plot.
%
% INPUT:
%   cfg                Configuration from reproduction_config(); uses
%                        rpo_to_deep_space.shared and diagnostics.
%   phase_database     Phase-aware RPO-Moon table from rpo2ds.load_database.
%   altitude_database  Historical dense or generated table for patch searches.
%   targets            Venus/Mars input struct from rpo2ds.load_targets.
%
% OUTPUT:
%   result             Struct with data.encounter_map, targets, target_transfers,
%                        earth_departure, altitude_sensitivity, summary metrics,
%                        and database-role provenance. Field suffixes identify
%                        km, km/s, days, and rad; no transfer NLP is solved.

settings = cfg.rpo_to_deep_space;
result = struct();
overall_timer = tic;

% Encounter map: retain a deterministic subset for plotting, but calculate
% every range and count from the complete phase-aware table.
map = struct();
map.time_of_flight_days = phase_database.t_event_days(:);
map.vinf_moon_km_s = phase_database.vInfMagMoon(:);
map.moon_soi_entries = phase_database.n_hits(:);
map.source_rows = local_plot_indices(map, settings.diagnostics.encounter_sample_limit);
result.data.encounter_map = map;
result.metrics.encounter_map = struct( ...
    'row_count', height(phase_database), ...
    'time_of_flight_range_days', [min(map.time_of_flight_days), max(map.time_of_flight_days)], ...
    'vinf_moon_range_km_s', [min(map.vinf_moon_km_s), max(map.vinf_moon_km_s)], ...
    'soi_entry_categories', unique(map.moon_soi_entries(:), 'sorted').');

result.data.targets = local_target_summary(targets);
result.metrics.targets = local_target_metrics(targets);
result.data.target_transfers = local_target_transfers(targets);

% Build one explanatory Earth-perigee construction from the Venus target.
result.data.earth_departure = local_build_earth_departure(settings, targets.venus);

if ismember('orbitPhase', altitude_database.Properties.VariableNames)
    altitude_database_role = "phase_aware_generated";
else
    altitude_database_role = "legacy_dense_reference";
end
% Compare two escape directions at fixed speed to show why the perigee
% altitude must be searched when constructing a patched initial guess.
comparison = local_compare_escape_directions( ...
    settings.diagnostics, altitude_database, altitude_database_role);
result.data.altitude_sensitivity = struct( ...
    'altitude_km', comparison.altitude_km, ...
    'direction_comparison', rmfield(comparison, 'altitude_km'), ...
    'database_role', altitude_database_role);
result.metrics.altitude_direction_comparison = struct( ...
    'database_role', altitude_database_role, ...
    'best_delta_v_range_km_s', local_finite_range( ...
        comparison.best_delta_v_km_s), ...
    'candidate_count_range', local_finite_range( ...
        comparison.n_candidates_in_rpo_vinf_range));

result.provenance.database_roles = struct( ...
    'encounter_map', "phase_aware", ...
    'altitude_sensitivity', altitude_database_role);
result.provenance.repository_root_used_for_inputs = false;
result.provenance.elapsed_seconds = toc(overall_timer);
end

function indices = local_plot_indices(map, sample_limit)
%LOCAL_PLOT_INDICES Select deterministic rows for encounter-map rendering.
%
% PROTOTYPE:
%   indices = local_plot_indices(map, sample_limit)
%
% DESCRIPTION:
%   Bin flight time and lunar speed, then include representatives of the
%   Moon-SOI hit categories before applying the final sample cap. Full-table
%   counts and ranges are computed elsewhere, independently of this subset.
%
% INPUT:
%   map           Struct with time_of_flight_days, vinf_moon_km_s,
%                 and moon_soi_entries arrays of equal length.
%   sample_limit  Maximum number of rows retained for rendering.
%
% OUTPUT:
%   indices       Column of source-row indices, in deterministic order.

n_rows = numel(map.time_of_flight_days);
if n_rows <= sample_limit
    indices = (1:n_rows).';
    return
end

% Deterministic binning preserves every rare n_hits category and coverage in
% the two plotted continuous quantities without storing the full database.
n_bins_t = 100;
n_bins_v = 20;
t = map.time_of_flight_days;
v = map.vinf_moon_km_s;
t_bin = min(n_bins_t, max(1, ceil((t - min(t)) / max(eps, range(t)) * n_bins_t)));
v_bin = min(n_bins_v, max(1, ceil((v - min(v)) / max(eps, range(v)) * n_bins_v)));
key = t_bin + n_bins_t * (v_bin - 1);
[~, first] = unique(key, 'stable');
indices = first(:);
for category = unique(map.moon_soi_entries(:)).'
    indices(end + 1, 1) = find(map.moon_soi_entries == category, 1, 'first'); %#ok<AGROW>
end
indices = unique(indices, 'stable');
if numel(indices) > sample_limit
    indices = indices(round(linspace(1, numel(indices), sample_limit))).';
end
end

function target_summary = local_target_summary(targets)
%LOCAL_TARGET_SUMMARY Keep the heliocentric target states needed by the diagnostic figures.

for name = ["venus", "mars"]
    target = targets.(name);
    target_summary.(name) = struct( ...
        'body_name', target.body_name, ...
        'departure_position_km', target.earth_departure_position_km(:), ...
        'arrival_position_km', target.planet_arrival_position_km(:), ...
        'departure_velocity_km_s', target.lambert_departure_velocity_km_s(:));
end
end

function metrics = local_target_metrics(targets)
%LOCAL_TARGET_METRICS Collect accepted epochs, flight times, and Earth-relative v-infinity vectors.

for name = ["venus", "mars"]
    target = targets.(name);
    metrics.(name) = struct( ...
        'departure_mjd2000_days', target.departure_mjd2000_days, ...
        'arrival_mjd2000_days', target.arrival_mjd2000_days, ...
        'time_of_flight_days', target.time_of_flight_days, ...
        'vinf_earth_magnitude_km_s', target.vinf_earth_magnitude_km_s, ...
        'vinf_earth_eci_km_s', target.vinf_earth_eci_km_s(:), ...
        'vinf_earth_rth_km_s', target.vinf_earth_rth_km_s(:));
end
end

function transfers = local_target_transfers(targets)
%LOCAL_TARGET_TRANSFERS Propagate the stored heliocentric Lambert departures for display.

mu_sun = targets.sun_gravitational_parameter_km3_s2;
for name = ["venus", "mars"]
    target = targets.(name);
    [~, arrival_velocity] = rpo2ds.solve_lambert( ...
        target.earth_departure_position_km, target.planet_arrival_position_km, ...
        target.time_of_flight_days, mu_sun);
    transfers.(name) = struct( ...
        'earth_position_km', target.earth_departure_position_km(:), ...
        'planet_position_km', target.planet_arrival_position_km(:), ...
        'departure_velocity_km_s', target.lambert_departure_velocity_km_s(:), ...
        'arrival_velocity_km_s', arrival_velocity(:), ...
        'positions_km', local_two_body_arc(target.earth_departure_position_km, ...
            target.lambert_departure_velocity_km_s, target.time_of_flight_days * 86400, mu_sun));
end
end

function earth_departure = local_build_earth_departure(settings, target)
%LOCAL_BUILD_EARTH_DEPARTURE Construct the explanatory Earth-perigee departure geometry.
%
% PROTOTYPE:
%   earth_departure = local_build_earth_departure(settings, target)
%
% DESCRIPTION:
%   Sample the perigee locus for the accepted Venus target, retain three
%   illustrative outgoing hyperbolae, and propagate representative Moon-Earth
%   Lambert legs selected near the fixed 1, 2, 4, 8, and 15 day targets.
%
% INPUT:
%   settings  Run02 study settings; shared contains the physical scales/body GM.
%   target    Accepted Venus target with outgoing Earth ECI v-infinity [km/s].
%
% OUTPUT:
%   earth_departure  Struct with perigee locus/velocities, hyperbola samples,
%                      and Lambert-arc cells in inertial km/km/s and days.

shared = settings.shared;
mu_earth = shared.earth_gravitational_parameter_km3_s2;
perigee_radius = shared.earth_radius_km + shared.minimum_earth_altitude_km;
[perigee_positions, perigee_velocities] = rp_circle_from_vinfplus( ...
    target.vinf_earth_eci_km_s, perigee_radius, mu_earth, 101);

earth_departure.perigee_locus_km = perigee_positions;
earth_departure.perigee_velocity_km_s = perigee_velocities;
earth_departure.outgoing_hyperbolae = local_hyperbolae(perigee_positions, ...
    perigee_velocities, mu_earth);

% Select five deterministic transfer-time quantiles from a small instance of
% the same Moon--Earth patched-conic database used by the transfer method.
moon_earth = interplanetary_MoonEarth_patch_opt( ...
    target.vinf_earth_eci_km_s, 11, shared.minimum_earth_altitude_km, ...
    2 * pi, 0, 20);
representative_times_days = [1, 2, 4, 8, 15];
arc_indices = zeros(size(representative_times_days));
for index = 1:numel(representative_times_days)
    [~, arc_indices(index)] = min(abs( ...
        moon_earth.ToF_days - representative_times_days(index)));
end
earth_departure.lambert_arc_indices = arc_indices(:);
earth_departure.lambert_arcs = cell(numel(arc_indices), 1);
for index = 1:numel(arc_indices)
    row = arc_indices(index);
    moon_position = moon_earth.R1(row, :);
    perigee_position = moon_earth.R2(row, :);
    departure_velocity = moon_earth.VI(row, :);
    arrival_velocity = moon_earth.VF(row, :);
    tof_seconds = moon_earth.ToF(row);
    arc = struct();
    arc.moon_position_km = moon_position;
    arc.perigee_position_km = perigee_position;
    arc.time_of_flight_days = tof_seconds / 86400;
    arc.departure_velocity_km_s = departure_velocity;
    arc.arrival_velocity_km_s = arrival_velocity;
    arc.perigee_velocity_km_s = moon_earth.Vp(row, :);
    arc.positions_km = local_two_body_arc(moon_position, departure_velocity, ...
        tof_seconds, mu_earth);
    earth_departure.lambert_arcs{index} = arc;
end
end

function arcs = local_hyperbolae(perigee_positions, perigee_velocities, mu)
%LOCAL_HYPERBOLAE Sample three outgoing hyperbolae and retain their near-Earth portions.

indices = round(linspace(1, size(perigee_positions, 1), 3));
arcs = cell(numel(indices), 1);
for index = 1:numel(indices)
    row = indices(index);
    positions = local_two_body_arc(perigee_positions(row, :), ...
        perigee_velocities(row, :), 3 * 86400, mu);
    % Keep the explanatory panel at Earth scale; the full escape arc is
    % represented by the heliocentric target figure, not this construction.
    near_earth = vecnorm(positions, 2, 2) <= 50000;
    if nnz(near_earth) < 2
        near_earth(1:min(10, size(positions, 1))) = true;
    end
    arcs{index} = struct('index', row, 'positions_km', positions(near_earth, :));
end
end

function positions_km = local_two_body_arc(position_km, velocity_km_s, duration_seconds, mu)
%LOCAL_TWO_BODY_ARC Return 200 sampled positions [km] over the requested duration [s].

sample_times = linspace(0, duration_seconds, 200).';
state0 = [position_km(:); velocity_km_s(:)];
solver_options = odeset('RelTol', 1e-10, 'AbsTol', 1e-11);
[~, states] = ode113(@(time, state) local_two_body_rhs(time, state, mu), ...
    sample_times, state0, solver_options);
positions_km = states(:, 1:3);
end

function derivative = local_two_body_rhs(~, state, mu)
%LOCAL_TWO_BODY_RHS Evaluate [v; -mu*r/|r|^3] for Cartesian km, km/s, and seconds.

position = state(1:3);
velocity = state(4:6);
derivative = [velocity; -mu * position / norm(position)^3];
end

function comparison = local_compare_escape_directions(settings, altitude_database, altitude_database_role)
%LOCAL_COMPARE_ESCAPE_DIRECTIONS Compare patched costs over the frozen perigee-altitude grid.
%
% PROTOTYPE:
%   comparison = local_compare_escape_directions(settings, altitude_database, ...
%       altitude_database_role)
%
% DESCRIPTION:
%   Evaluate each altitude independently for both selected escape directions
%   at a common speed. Retain candidate counts even when no patch is selected;
%   use the same per-altitude calculation in serial and parallel loops.
%
% INPUT:
%   settings                Run02 diagnostics settings and sampling grid.
%   altitude_database       RPO-Moon table used for the patch search.
%   altitude_database_role  Source-role label retained in diagnostics.
%
% OUTPUT:
%   comparison              Altitude [km], speed [km/s], directions, best
%                             perigee costs [km/s], and candidate-count arrays.

altitudes = settings.altitude_samples(:);
speed = settings.direction_comparison_speed_km_s;
directions = [settings.outgoing_direction(:).'; settings.direction_comparison_vector];
best_delta_v = NaN(numel(altitudes), size(directions, 1));
n_candidates = NaN(size(best_delta_v));
fprintf('Escape-direction comparison uses %s database at %.1f km/s.\n', ...
    altitude_database_role, speed);
timer = tic;
if rpo.parallel_available(settings.use_parallel)
    fprintf('Parallel Computing Toolbox detected; the comparison will use parallel execution.\n');
    completed_count = 0;
    report_stride = max(1, ceil(numel(altitudes) / 10));
    progress_queue = parallel.pool.DataQueue;
    afterEach(progress_queue, @report_parallel_progress);
    parfor altitude_index = 1:numel(altitudes)
        [best_delta_v(altitude_index, :), n_candidates(altitude_index, :)] = ...
            local_altitude_row(altitudes(altitude_index), speed, ...
            directions, settings, altitude_database);
        send(progress_queue, 1);
    end
else
    fprintf(['Parallel Computing Toolbox is not available. The comparison ', ...
        'will run serially and take longer; install the toolbox for faster execution.\n']);
    for altitude_index = 1:numel(altitudes)
        [best_delta_v(altitude_index, :), n_candidates(altitude_index, :)] = ...
            local_altitude_row(altitudes(altitude_index), speed, ...
            directions, settings, altitude_database);
        fprintf('  Escape-direction progress: %d/%d (%.1f s elapsed)\n', ...
            altitude_index, numel(altitudes), toc(timer));
    end
end
comparison = struct('altitude_km', altitudes, ...
    'speed_km_s', speed, 'directions', directions, ...
    'best_delta_v_km_s', best_delta_v, ...
    'n_candidates_in_rpo_vinf_range', n_candidates, ...
    'database_role', altitude_database_role);

    function report_parallel_progress(~)
    %REPORT_PARALLEL_PROGRESS Count completed workers and print campaign progress.

        completed_count = completed_count + 1;
        if mod(completed_count, report_stride) == 0 || completed_count == numel(altitudes)
            fprintf('  Escape-direction progress: %d/%d (%.1f s elapsed)\n', ...
                completed_count, numel(altitudes), toc(timer));
        end
    end
end

function [best_delta_v_row, n_candidates_row] = local_altitude_row(altitude_km, speed, directions, settings, altitude_database)
%LOCAL_ALTITUDE_ROW Search Moon-Earth/RPO-Moon patches for both directions at one altitude [km].

best_delta_v_row = NaN(1, size(directions, 1));
n_candidates_row = NaN(size(best_delta_v_row));
for direction_index = 1:size(directions, 1)
    vinf = speed * directions(direction_index, :);
    moon_earth_database = interplanetary_MoonEarth_patch_opt( ...
        vinf, settings.perigee_samples, altitude_km, 2 * pi, 0, 100);
    [solution, diagnostics] = MoonEarth_RPOMoon_cellVersion_opt( ...
        moon_earth_database, altitude_database, settings.best_solution_count, ...
        settings.lunar_vinf_tolerance_km_s, [], struct('mode', 'optimization_only'));
    if RPO2DS_has_patch_solution(solution)
        best_delta_v_row(direction_index) = solution.deltavpMag(1);
    end
    n_candidates_row(direction_index) = diagnostics.vInfMoon_within_constraints;
end
end

function values = local_finite_range(array)
%LOCAL_FINITE_RANGE Return finite extrema, or [NaN, NaN] when no finite values remain.

finite_values = array(isfinite(array));
if isempty(finite_values)
    values = [NaN, NaN];
else
    values = [min(finite_values), max(finite_values)];
end
end
