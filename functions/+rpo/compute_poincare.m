function result = compute_poincare(cfg, a1_result)
%COMPUTE_POINCARE Sample the phase-space structure around the reference RPO.
%
% PROTOTYPE:
%   result = rpo.compute_poincare(cfg, a1_result)
%
% DESCRIPTION:
%   Draw seeded initial conditions uniformly within an (x, vx) ellipse at
%   y = 0 and the reference Jacobi constant. Retain only complete chains
%   of positive y = 0 crossings. rng is reset to the configured seed.
%
% INPUT:
%   cfg         Configuration from reproduction_config(); uses constants
%               and reference_rpo.poincare settings.
%   a1_result   A1 result supplying data.rpo.state [4 x 1], ordered
%               [x; y; vx; vy] in the Earth-Moon rotating frame.
%
% OUTPUT:
%   result      Struct with initial_conditions [N x 4] and retained crossing
%               arrays [M x K] in data; M is the number of valid chains and
%               K the requested crossing count. x is in DU, vx in VU, and
%               crossing times in TU. Metrics count valid/rejected chains.
%               Configuration and provenance are attached separately.

settings = cfg.reference_rpo.poincare;
mu = cfg.constants.mu_earth_moon;
reference_state = a1_result.data.rpo.state;
jacobi_constant = rpo.jacobi_constant(reference_state', mu);

% Draw uniformly in ellipse area, then recover positive vy from the Jacobi integral.
rng(settings.random_seed, 'twister');
angles = 2 * pi * rand(settings.initial_condition_count, 1);
radii = sqrt(rand(settings.initial_condition_count, 1));
x0 = reference_state(1) + settings.ellipse_x_radius_du ...
    .* radii .* cos(angles);
vx0 = reference_state(3) + settings.ellipse_vx_radius_vu ...
    .* radii .* sin(angles);
potential = local_potential(x0, zeros(size(x0)), mu);
vy0 = sqrt(2 * potential - jacobi_constant - vx0 .^ 2);
initial_conditions = [x0, zeros(size(x0)), vx0, vy0];

% Allocate complete crossing chains and choose the execution mode.
n_initial = size(initial_conditions, 1);
n_crossings = settings.crossing_count;
crossing_x = NaN(n_initial, n_crossings);
crossing_vx = NaN(n_initial, n_crossings);
crossing_times = NaN(n_initial, n_crossings);
valid_chains = false(n_initial, 1);
use_parallel = rpo.parallel_available(settings.use_parallel);

% Propagate every initial condition through the requested positive crossings.
fprintf('Poincare section: %d initial conditions, %d crossings each, %s execution.\n', ...
    n_initial, n_crossings, local_execution_name(use_parallel));
if use_parallel
    parfor index = 1:n_initial
        [crossing_x(index, :), crossing_vx(index, :), ...
            crossing_times(index, :), valid_chains(index)] = ...
            local_chain(initial_conditions(index, :), settings, mu);
    end
else
    for index = 1:n_initial
        [crossing_x(index, :), crossing_vx(index, :), ...
            crossing_times(index, :), valid_chains(index)] = ...
            local_chain(initial_conditions(index, :), settings, mu);
    end
end

% Keep complete chains for plotting and retain the rejection mask for diagnostics.
result.data.reference_state = reference_state;
result.data.jacobi_constant = jacobi_constant;
result.data.initial_conditions = initial_conditions;
result.data.crossing_x = crossing_x(valid_chains, :);
result.data.crossing_vx = crossing_vx(valid_chains, :);
result.data.crossing_times_tu = crossing_times(valid_chains, :);
result.data.valid_chain_mask = valid_chains;
result.metrics.requested_initial_condition_count = n_initial;
result.metrics.valid_chain_count = nnz(valid_chains);
result.metrics.rejected_chain_count = nnz(~valid_chains);
result.metrics.crossings_per_valid_chain = n_crossings;
end

function [x, vx, times, valid] = local_chain(initial_condition, settings, mu)
%LOCAL_CHAIN Return one complete chain of positive y = 0 crossings.
%
% initial_condition is [x, y, vx, vy] in DU/VU. settings supplies the
% crossing count and integration limits. Outputs x, vx, and times are rows
% in DU, VU, and cumulative TU; valid requires the full chain.

[crossings, times, valid] = rpo.propagate_crossings( ...
    initial_condition, settings.crossing_count, ...
    settings.crossing_time_limit_tu, mu, false, settings.event_tolerance, ...
    1, true);
x = NaN(1, settings.crossing_count);
vx = NaN(1, settings.crossing_count);
if valid
    for crossing = 1:settings.crossing_count
        x(crossing) = crossings{crossing}.state(1);
        vx(crossing) = crossings{crossing}.state(3);
    end
end
times = times';
end

function potential = local_potential(x, y, mu)
%LOCAL_POTENTIAL Evaluate normalized CR3BP effective potential U(x, y).
%
% x and y are equal-size position arrays [DU]; mu is the mass ratio [-].

r1 = hypot(x + mu, y);
r2 = hypot(x - 1 + mu, y);
potential = 0.5 * (x .^ 2 + y .^ 2) + (1 - mu) ./ r1 + mu ./ r2;
end

function name = local_execution_name(use_parallel)
%LOCAL_EXECUTION_NAME Return the serial or parallel label for progress output.

if use_parallel
    name = 'parallel';
else
    name = 'serial';
end
end
