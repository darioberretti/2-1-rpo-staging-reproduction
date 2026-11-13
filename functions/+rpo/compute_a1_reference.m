function result = compute_a1_reference(cfg)
%COMPUTE_A1_REFERENCE Reconstruct the A1 family and corrected 2:1 RPO.
%
% PROTOTYPE:
%   result = rpo.compute_a1_reference(cfg)
%
% DESCRIPTION:
%   Sweep and refine symmetric-orbit seeds, continue the three A1 branches,
%   and correct the reference orbit to the prescribed synodic resonance.
%   States use the Earth-Moon rotating frame and normalized DU, TU, and VU.
%
% INPUT:
%   cfg      Configuration from reproduction_config(); uses constants
%            and reference_rpo.a1 / reference_rpo.rpo settings.
%
% OUTPUT:
%   result   Struct with data, metrics, and convergence fields.
%            data.family is [N x 5]: [x0, vy0, half-period, final x, k].
%            Positions are in DU, velocities in VU, and times in TU.
%            data.rpo contains state [4 x 1] = [x; y; vx; vy], period,
%            and sampled trajectory. Metrics include closure and Moon altitude.
%            Configuration and provenance are attached by rpo.finalize_result.

% Read the frozen sweep settings and choose serial or parallel execution.
a1 = cfg.reference_rpo.a1;
mu = cfg.constants.mu_earth_moon;
use_parallel = rpo.parallel_available(a1.use_parallel);

% Sweep initial conditions and record successive y = 0 crossings.
fprintf('A1 sweep: %d x %d points using %s execution.\n', ...
    numel(a1.x0_values), numel(a1.vy0_values), ...
    local_execution_name(use_parallel));
[crossing_vx, crossing_times, sweep_rejections] = local_sweep(a1, mu, use_parallel);

% Find sign-change seeds, correct their symmetry, and merge duplicates.
seeds = rpo.collect_a1_seeds(a1.x0_values, a1.vy0_values, ...
    crossing_times, crossing_vx, a1.seed_velocity_merge_tolerance);
fprintf('Refining %d sign-change seeds.\n', numel(seeds));
[refined_seeds, refinement_rejections] = local_refine(seeds, a1, mu, use_parallel);
unique_seeds = local_merge_refined(refined_seeds, a1.seed_velocity_merge_tolerance);

% Continue each physical branch from its frozen starting neighborhood.
k1_start = local_select_seed(unique_seeds, a1.branches.k1);
k2_start = local_select_seed(unique_seeds, a1.branches.k2);
k3_start = local_select_seed(unique_seeds, a1.branches.k3);

k1 = local_continue(a1.branches.k1.x0_values, k1_start.state(4), 1, a1, mu);
k1 = local_retain(k1, a1.branches.k1.retained_x0_range);

k2_forward = local_continue(a1.branches.k2.forward_x0_values, ...
    k2_start.state(4), 2, a1, mu);
k2_backward = local_continue(a1.branches.k2.backward_x0_values, ...
    k2_start.state(4), 2, a1, mu);
k2 = local_retain([flipud(k2_backward); k2_forward], ...
    a1.branches.k2.retained_x0_range);

k3_forward = local_continue(a1.branches.k3.forward_x0_values, ...
    k3_start.state(4), 3, a1, mu);
k3_backward = local_continue(a1.branches.k3.backward_x0_values, ...
    k3_start.state(4), 3, a1, mu);
k3 = flipud(local_retain([flipud(k3_backward); k3_forward], ...
    a1.branches.k3.retained_x0_range));

% Assemble the family, correct the resonant RPO, and sample plot trajectories.
family = [k1; k2; k3];
rpo_orbit = local_correct_rpo(family, cfg.reference_rpo.rpo, mu);
[sample_rows, sample_trajectories] = local_sample_trajectories( ...
    family, a1.family_sample_indices, mu);

% Propagate the reference orbit to measure closure and half-period Moon distance.
period = 2 * rpo_orbit.t_half;
[final_state, ~, ~, trajectory, trajectory_times] = propagate( ...
    rpo_orbit.state, 0, period, mu, '2D', 'CR3BP', 'false', false, ...
    [], [], cfg.reference_rpo.rpo.integration_tolerance, ...
    cfg.reference_rpo.rpo.integration_tolerance);
moon_position = [1 - mu, 0];
half_state = propagate(rpo_orbit.state, 0, rpo_orbit.t_half, mu, ...
    '2D', 'CR3BP', 'false', false, [], [], ...
    cfg.reference_rpo.rpo.integration_tolerance, ...
    cfg.reference_rpo.rpo.integration_tolerance);
minimum_moon_distance_km = norm(half_state(1:2)' - moon_position) ...
    * cfg.constants.distance_unit_km;

% Keep the sweep, branch, and trajectory data needed for replay.
result.data.sweep.x0_values = a1.x0_values;
result.data.sweep.vy0_values = a1.vy0_values;
result.data.sweep.crossing_vx = crossing_vx;
result.data.sweep.crossing_times = crossing_times;
result.data.seeds = seeds;
result.data.refined_seeds = refined_seeds;
result.data.unique_seeds = unique_seeds;
result.data.branches.k1 = k1;
result.data.branches.k2 = k2;
result.data.branches.k3 = k3;
result.data.family = family;
result.data.rpo.state = rpo_orbit.state;
result.data.rpo.period_tu = period;
result.data.rpo.trajectory = trajectory;
result.data.rpo.trajectory_times_tu = trajectory_times;
result.data.figure_data.family_sample_indices = a1.family_sample_indices;
result.data.figure_data.family_sample_rows = sample_rows;
result.data.figure_data.family_sample_trajectories = sample_trajectories;
result.data.figure_data.rpo_neighbor_count = a1.rpo_neighbor_count;

% Summarize branch ranges and the reference orbit in normalized/physical units.
result.metrics.branch_ranges.k1 = local_branch_metrics(k1, mu);
result.metrics.branch_ranges.k2 = local_branch_metrics(k2, mu);
result.metrics.branch_ranges.k3 = local_branch_metrics(k3, mu);
result.metrics.rpo.state = rpo_orbit.state;
result.metrics.rpo.period_tu = period;
result.metrics.rpo.period_days = period * cfg.constants.time_unit_days;
result.metrics.rpo.half_period_tu = rpo_orbit.t_half;
result.metrics.rpo.half_period_days = rpo_orbit.t_half ...
    * cfg.constants.time_unit_days;
result.metrics.rpo.jacobi_constant = rpo.jacobi_constant(rpo_orbit.state', mu);
result.metrics.rpo.minimum_moon_distance_km = minimum_moon_distance_km;
result.metrics.rpo.minimum_moon_altitude_km = minimum_moon_distance_km ...
    - cfg.constants.moon_radius_km;
result.metrics.rpo.closure_residual = norm(final_state - rpo_orbit.state);

% Record rejected integrations and correction convergence diagnostics.
result.convergence.sweep_rejections = sweep_rejections;
result.convergence.refinement_rejections = refinement_rejections;
result.convergence.refined_seed_count = nnz([refined_seeds.valid]);
result.convergence.unique_seed_count = numel(unique_seeds);
result.convergence.rpo_iterations = rpo_orbit.iterations;
result.convergence.rpo_residual = rpo_orbit.residual;
fprintf('Rejected candidate integrations: %d sweep, %d refinement.\n', ...
    sweep_rejections, refinement_rejections);
end

function [velocities, times, rejection_count] = local_sweep(cfg, mu, use_parallel)
%LOCAL_SWEEP Evaluate crossing vx and time over the coarse A1 grid.
%
% PROTOTYPE:
%   [velocities, times, rejection_count] = local_sweep(cfg, mu, use_parallel)
%
% DESCRIPTION:
%   Process one x0 column at a time and optionally parallelize the vy0 rows.
%   Invalid crossing chains remain NaN; progress is reported by x0 column.
%
% INPUT:
%   cfg           Frozen reference_rpo.a1 settings.
%   mu            Secondary-to-total mass ratio [-].
%   use_parallel  Logical selecting the parallel or serial inner loop.
%
% OUTPUT:
%   velocities       [Nv x Nx x K] Crossing vx [VU].
%   times            [Nv x Nx x K] Crossing times [TU].
%   rejection_count  Number of grid points with an invalid crossing chain.
%                      Nv, Nx, and K count velocities, positions, and crossings.

n_x = numel(cfg.x0_values);
n_vy = numel(cfg.vy0_values);
n_crossings = numel(cfg.crossing_numbers);
vy0_values = cfg.vy0_values;
velocities = NaN(n_vy, n_x, n_crossings);
times = NaN(n_vy, n_x, n_crossings);
rejection_count = 0;
start = tic;

% Process one x0 column at a time; parallel workers span the vy0 grid.
for x_index = 1:n_x
    column_velocities = NaN(n_vy, n_crossings);
    column_times = NaN(n_vy, n_crossings);
    column_rejected = false(n_vy, 1);
    x0 = cfg.x0_values(x_index);
    if use_parallel
        parfor vy_index = 1:n_vy
            [column_velocities(vy_index, :), column_times(vy_index, :), ...
                column_rejected(vy_index)] = local_sweep_point( ...
                x0, vy0_values(vy_index), cfg, mu);
        end
    else
        for vy_index = 1:n_vy
            [column_velocities(vy_index, :), column_times(vy_index, :), ...
                column_rejected(vy_index)] = local_sweep_point( ...
                x0, vy0_values(vy_index), cfg, mu);
        end
    end

    % Store the completed column and report coarse-grained progress.
    velocities(:, x_index, :) = column_velocities;
    times(:, x_index, :) = column_times;
    rejection_count = rejection_count + nnz(column_rejected);
    if mod(x_index, cfg.progress_interval) == 0 || x_index == n_x
        local_progress('A1 sweep', x_index, n_x, start);
    end
end
end

function [velocities, times, rejected] = local_sweep_point(x0, vy0, cfg, mu)
%LOCAL_SWEEP_POINT Evaluate one [x0, 0, 0, vy0] crossing chain.
%
% Returns row vectors of crossing vx [VU] and times [TU]; rejected is true
% if the complete chain cannot be accepted.

n_crossings = numel(cfg.crossing_numbers);
velocities = NaN(1, n_crossings);
times = NaN(1, n_crossings);

% Retain all requested crossings only when the complete chain is valid.
[crossings, crossing_times, valid] = rpo.propagate_crossings( ...
    [x0; 0; 0; vy0], n_crossings, cfg.coarse_time_limit_tu, mu, false, ...
    cfg.event_tolerance);
if valid
    for crossing = 1:n_crossings
        velocities(crossing) = crossings{crossing}.state(3);
        times(crossing) = crossing_times(crossing);
    end
end
rejected = ~valid;
end

function [refined, rejection_count] = local_refine(seeds, cfg, mu, use_parallel)
%LOCAL_REFINE Differentially correct the coarse sign-change seeds.
%
% PROTOTYPE:
%   [refined, rejection_count] = local_refine(seeds, cfg, mu, use_parallel)
%
% DESCRIPTION:
%   Correct seeds in fixed-size chunks, preserving their order and crossing
%   numbers. The serial and parallel paths use the same correction routine.
%
% INPUT:
%   seeds         Struct array from rpo.collect_a1_seeds.
%   cfg           Frozen reference_rpo.a1 correction settings.
%   mu            Secondary-to-total mass ratio [-].
%   use_parallel  Logical selecting the parallel or serial loop.
%
% OUTPUT:
%   refined          Struct array of rpo.correct_symmetric_orbit outputs,
%                      each with its original crossing_number attached.
%   rejection_count  Number of corrected candidates with valid = false.

empty = struct('valid', false, 'state', NaN(4, 1), 't_half', NaN, ...
    'x_final', NaN, 'iterations', 0, 'residual', Inf, ...
    'crossing_number', 0);
refined = repmat(empty, numel(seeds), 1);
chunk_size = 250;
start = tic;

% Correct seeds in chunks so progress remains visible during parallel runs.
for first = 1:chunk_size:numel(seeds)
    last = min(first + chunk_size - 1, numel(seeds));
    if use_parallel
        parfor index = first:last
            orbit = rpo.correct_symmetric_orbit(seeds(index).x0, ...
                seeds(index).vy0, seeds(index).crossing_number, cfg, mu);
            orbit.crossing_number = seeds(index).crossing_number;
            refined(index) = orbit;
        end
    else
        for index = first:last
            orbit = rpo.correct_symmetric_orbit(seeds(index).x0, ...
                seeds(index).vy0, seeds(index).crossing_number, cfg, mu);
            orbit.crossing_number = seeds(index).crossing_number;
            refined(index) = orbit;
        end
    end
    local_progress('A1 refinement', last, numel(seeds), start);
end
rejection_count = nnz(~[refined.valid]);
end

function unique_seeds = local_merge_refined(refined, merge_tolerance)
%LOCAL_MERGE_REFINED Retain distinct valid corrections at each x0.
%
% PROTOTYPE:
%   unique_seeds = local_merge_refined(refined, merge_tolerance)
%
% DESCRIPTION:
%   Preserve the first accepted candidate for each crossing number and
%   initial x0 within 1e-12 DU, merging initial vy within the given tolerance.
%
% INPUT:
%   refined          Struct array from local_refine.
%   merge_tolerance  Maximum initial vy separation treated as duplicate [VU].
%
% OUTPUT:
%   unique_seeds     Valid, de-duplicated orbit structs in input traversal order.


% Discard failed corrections before comparing states at the same crossing number.
refined = refined([refined.valid]);
unique_seeds = refined([]);
for index = 1:numel(refined)
    candidate = refined(index);
    duplicate = false;
    for accepted = 1:numel(unique_seeds)
        duplicate = unique_seeds(accepted).crossing_number ...
            == candidate.crossing_number ...
            && abs(unique_seeds(accepted).state(1) - candidate.state(1)) < 1e-12 ...
            && abs(unique_seeds(accepted).state(4) - candidate.state(4)) ...
            <= merge_tolerance;
        if duplicate
            break
        end
    end
    if ~duplicate
        unique_seeds(end + 1) = candidate; %#ok<AGROW>
    end
end
end

function seed = local_select_seed(seeds, branch)
%LOCAL_SELECT_SEED Select the seed nearest a frozen physical branch start.
%
% PROTOTYPE:
%   seed = local_select_seed(seeds, branch)
%
% DESCRIPTION:
%   Restrict candidates to the branch crossing number and rank normalized
%   (x0, vy0) distance using the frozen 0.01 scales. Fail if no candidate
%   lies within the implemented distance threshold.
%
% INPUT:
%   seeds   Valid corrected seed structs from local_merge_refined.
%   branch  A1 branch settings containing crossing_number and start [x0, vy0].
%
% OUTPUT:
%   seed    Selected orbit struct; x0 is in DU and vy0 in VU.

matching = find([seeds.crossing_number] == branch.crossing_number);
if isempty(matching)
    error('rpo:A1BranchSeedMissing', ...
        'No refined seed is available for A1 branch k=%d.', branch.crossing_number);
end

% Rank candidates by distance to the frozen physical branch start.
states = reshape([seeds(matching).state], 4, [])';
distance = hypot((states(:, 1) - branch.start(1)) / 0.01, ...
    (states(:, 4) - branch.start(2)) / 0.01);
[minimum_distance, local_index] = min(distance);
if minimum_distance > 10
    error('rpo:A1BranchSeedMissing', ...
        'No refined seed is sufficiently close to the physical k=%d start.', ...
        branch.crossing_number);
end
seed = seeds(matching(local_index));
end

function branch = local_continue(x0_values, initial_vy, crossing_number, cfg, mu)
%LOCAL_CONTINUE Continue symmetric orbits along a prescribed x0 grid.
%
% PROTOTYPE:
%   branch = local_continue(x0_values, initial_vy, crossing_number, cfg, mu)
%
% DESCRIPTION:
%   Use the most recent valid vy0 as the next initial guess. Failed grid
%   points remain NaN and do not replace the continuation guess.
%
% INPUT:
%   x0_values        Vector of initial x coordinates [DU], in traversal order.
%   initial_vy       Starting y-velocity guess [VU].
%   crossing_number  Selected y = 0 crossing number.
%   cfg              Frozen reference_rpo.a1 correction settings.
%   mu               Secondary-to-total mass ratio [-].
%
% OUTPUT:
%   branch           [N x 5] Rows [x0, vy0, half-period, final x, k], with
%                      positions in DU, velocity in VU, and time in TU.

branch = NaN(numel(x0_values), 5);
vy = initial_vy;

% Use each successful correction as the velocity guess for the next x0.
for index = 1:numel(x0_values)
    orbit = rpo.correct_symmetric_orbit(x0_values(index), vy, ...
        crossing_number, cfg, mu);
    if orbit.valid
        branch(index, :) = [orbit.state(1), orbit.state(4), orbit.t_half, ...
            orbit.x_final, crossing_number];
        vy = orbit.state(4);
    end
end
end

function branch = local_retain(branch, x0_range)
%LOCAL_RETAIN Keep finite branch rows within the frozen x0 range.
%
% branch uses [x0, vy0, half-period, final x, k]; x0_range is in DU.

valid = all(isfinite(branch), 2) & branch(:, 1) >= x0_range(1) - 1e-12 ...
    & branch(:, 1) <= x0_range(2) + 1e-12;
branch = branch(valid, :);
end

function orbit = local_correct_rpo(family, cfg, mu)
%LOCAL_CORRECT_RPO Correct the RPO at the prescribed fixed half-period.
%
% PROTOTYPE:
%   orbit = local_correct_rpo(family, cfg, mu)
%
% DESCRIPTION:
%   Choose the family member nearest the resonant half-period, then correct
%   x0 and vy0 using the fixed-time STM so y and vx vanish at that time.
%   Failure to meet both residual and step tolerances raises an error.
%
% INPUT:
%   family  [N x 5] Rows [x0, vy0, half-period, final x, k] in DU, VU, TU.
%   cfg     Frozen reference_rpo.rpo resonance and correction settings.
%   mu      Secondary-to-total mass ratio [-].
%
% OUTPUT:
%   orbit   Struct with state [4 x 1] = [x; y; vx; vy], t_half [TU],
%             iteration count, and maximum absolute symmetry residual.
%             The residual combines normalized y and vx components.


% Fix the resonant period and start from the nearest family member.
period = -pi / cfg.synodic_frequency_rad_tu;
t_half = period / 2;
[~, seed_index] = min(abs(family(:, 3) - t_half));
state = [family(seed_index, 1); 0; 0; family(seed_index, 2)];
previous_state = zeros(4, 1);
orbit = struct('state', state, 't_half', t_half, 'iterations', 0, ...
    'residual', Inf);

for iteration = 1:cfg.maximum_iterations
    [final, stm] = propagate(state, 0, t_half, mu, ...
        '2D', 'CR3BP', 'true', false, [], [], ...
        cfg.integration_tolerance, cfg.integration_tolerance);
    residual = -final(2:3);
    step = state - previous_state;
    previous_state = state;

    % Enforce y = vx = 0 at the fixed half-period by correcting x0 and vy0.
    correction = stm(2:3, [1, 4]) \ residual;
    state([1, 4]) = state([1, 4]) + correction;
    if all(abs(residual) <= cfg.residual_tolerance) ...
            && all(abs(step([1, 4])) <= cfg.step_tolerance)
        orbit.state = state;
        orbit.iterations = iteration;
        orbit.residual = max(abs(residual));
        return
    end
end
error('rpo:RPOCorrectionFailed', ...
    'The reference RPO correction did not converge in %d iterations.', ...
    cfg.maximum_iterations);
end

function metrics = local_branch_metrics(branch, mu)
%LOCAL_BRANCH_METRICS Summarize x0, Jacobi constant, and half-period ranges.
%
% branch rows are [x0, vy0, half-period, final x, k] in DU, VU, TU.

metrics.x0_range = [min(branch(:, 1)), max(branch(:, 1))];
states = [branch(:, 1), zeros(size(branch, 1), 1), ...
    zeros(size(branch, 1), 1), branch(:, 2)];
jacobi = rpo.jacobi_constant(states, mu);
metrics.jacobi_range = [min(jacobi), max(jacobi)];
metrics.half_period_range_tu = [min(branch(:, 3)), max(branch(:, 3))];
end

function local_progress(label, completed, total, start)
%LOCAL_PROGRESS Print completion percentage and an elapsed-time estimate.

elapsed = toc(start);
remaining = elapsed / completed * (total - completed);
fprintf('%s: %5.1f%%, elapsed %.1f min, remaining %.1f min.\n', ...
    label, 100 * completed / total, elapsed / 60, remaining / 60);
end

function name = local_execution_name(use_parallel)
%LOCAL_EXECUTION_NAME Return the serial or parallel label for progress output.

if use_parallel
    name = 'parallel';
else
    name = 'serial';
end
end

function [rows, trajectories] = local_sample_trajectories(family, family_indices, mu)
%LOCAL_SAMPLE_TRAJECTORIES Propagate selected family rows over one full period.
%
% family_indices selects rows [x0, vy0, half-period, final x, k].
% Returns the selected rows and trajectory cells in DU/VU; propagation
% uses the normalized rotating frame and the default integration tolerances.

rows = family(family_indices, :);
trajectories = cell(numel(family_indices), 1);
for index = 1:numel(family_indices)
    state = [rows(index, 1); 0; 0; rows(index, 2)];
    [~, ~, ~, trajectories{index}] = propagate(state, 0, ...
        2 * rows(index, 3), mu, '2D', 'CR3BP', 'false', false);
end
end
