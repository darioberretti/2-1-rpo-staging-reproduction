function result = compute_vertical(cfg, a1_result)
%COMPUTE_VERTICAL Propagate the frozen vertical-magnitude and phase studies.
%
% PROTOTYPE:
%   result = rpo.compute_vertical(cfg, a1_result)
%
% DESCRIPTION:
%   Embed the planar reference RPO in spatial coordinates, apply signed vz
%   impulses at the initial state and selected orbital phases, and propagate
%   each case over the frozen duration in the Earth-Moon rotating frame.
%
% INPUT:
%   cfg         Configuration from reproduction_config(); uses constants
%               and reference_rpo.vertical settings.
%   a1_result   A1 result supplying data.rpo.state [x; y; vx; vy] and
%               data.rpo.period_tu.
%
% OUTPUT:
%   result      Struct with magnitude and phase case data. Trajectory cells
%               contain rows [x, y, z, vx, vy, vz] in DU and VU; time cells
%               are supplied in TU and days. Phase fractions are [-].
%               Metrics contain peak absolute z in DU and km.
%               Configuration and provenance are attached separately.

% Embed the planar reference state in the spatial rotating-frame state order.
settings = cfg.reference_rpo.vertical;
constants = cfg.constants;
mu = constants.mu_earth_moon;
period = a1_result.data.rpo.period_tu;
reference_2d = a1_result.data.rpo.state;
reference_3d = [reference_2d(1:2); 0; reference_2d(3:4); 0];
phase_times = settings.phase_fractions * period;

% Recover the unperturbed states at the requested fractions of one RPO period.
phase_states = NaN(6, numel(phase_times));
for phase = 1:numel(phase_times)
    phase_states(:, phase) = propagate(reference_3d, 0, phase_times(phase), ...
        mu, '3D', 'CR3BP', 'false', false);
end

% Combine the signed-magnitude cases and the fixed-impulse phase cases.
n_magnitude_cases = numel(settings.impulses_vu);
n_phase_cases = numel(phase_times);
case_states = [repmat(reference_3d, 1, n_magnitude_cases), phase_states];
case_impulses = [settings.impulses_vu; ...
    repmat(settings.phase_impulse_vu, n_phase_cases, 1)];
n_cases = numel(case_impulses);
trajectories = cell(n_cases, 1);
times_tu = cell(n_cases, 1);
maximum_absolute_z_du = NaN(n_cases, 1);
duration_tu = settings.duration_days / constants.time_unit_days;
use_parallel = rpo.parallel_available(settings.use_parallel);

% Apply each vertical impulse and propagate over the common study duration.
fprintf(['Vertical response: %d signed-magnitude cases and %d phase cases ', ...
    'over %.0f days, %s execution.\n'], n_magnitude_cases, n_phase_cases, ...
    settings.duration_days, local_execution_name(use_parallel));
if use_parallel
    parfor case_index = 1:n_cases
        [trajectories{case_index}, times_tu{case_index}, ...
            maximum_absolute_z_du(case_index)] = local_propagate( ...
            case_states(:, case_index), case_impulses(case_index), ...
            duration_tu, mu);
    end
else
    for case_index = 1:n_cases
        [trajectories{case_index}, times_tu{case_index}, ...
            maximum_absolute_z_du(case_index)] = local_propagate( ...
            case_states(:, case_index), case_impulses(case_index), ...
            duration_tu, mu);
    end
end

% Separate the two case groups and store normalized and physical-unit outputs.
magnitude_indices = 1:n_magnitude_cases;
phase_indices = n_magnitude_cases + (1:n_phase_cases);
result.data.reference_state_3d = reference_3d;
result.data.magnitude.impulses_vu = settings.impulses_vu;
result.data.magnitude.trajectories = trajectories(magnitude_indices);
result.data.magnitude.times_tu = times_tu(magnitude_indices);
result.data.magnitude.times_days = local_days( ...
    times_tu(magnitude_indices), constants.time_unit_days);
result.data.phase.fractions = settings.phase_fractions;
result.data.phase.times_tu = phase_times;
result.data.phase.states = phase_states;
result.data.phase.impulse_vu = settings.phase_impulse_vu;
result.data.phase.trajectories = trajectories(phase_indices);
result.data.phase.propagation_times_tu = times_tu(phase_indices);
result.data.phase.propagation_times_days = local_days( ...
    times_tu(phase_indices), constants.time_unit_days);
result.metrics.maximum_absolute_z_by_impulse_du = ...
    maximum_absolute_z_du(magnitude_indices);
result.metrics.maximum_absolute_z_by_impulse_km = ...
    result.metrics.maximum_absolute_z_by_impulse_du ...
    * constants.distance_unit_km;
result.metrics.maximum_absolute_z_by_phase_du = ...
    maximum_absolute_z_du(phase_indices);
result.metrics.maximum_absolute_z_by_phase_km = ...
    result.metrics.maximum_absolute_z_by_phase_du ...
    * constants.distance_unit_km;
end

function days = local_days(times, time_unit_days)
%LOCAL_DAYS Convert each time-vector cell from TU to days.

days = cellfun(@(time) time * time_unit_days, times, ...
    'UniformOutput', false);
end

function [trajectory, time, maximum_absolute_z] = local_propagate( ...
        phase_state, impulse, duration, mu)
%LOCAL_PROPAGATE Apply a vz impulse and propagate one spatial case.
%
% phase_state is [x; y; z; vx; vy; vz] in DU/VU, impulse is in VU, and
% duration is in TU. Return trajectory rows in DU/VU, time in TU, and the
% maximum absolute sampled z in DU.

initial_state = phase_state;
initial_state(6) = initial_state(6) + impulse;
[~, ~, ~, trajectory, time] = propagate(initial_state, 0, duration, mu, ...
    '3D', 'CR3BP', 'false', false);
maximum_absolute_z = max(abs(trajectory(:, 3)));
end

function name = local_execution_name(use_parallel)
%LOCAL_EXECUTION_NAME Return the serial or parallel label for progress output.

if use_parallel
    name = 'parallel';
else
    name = 'serial';
end
end
