function result = compute_stability(cfg, a1_result)
%COMPUTE_STABILITY Compute planar and spatial A1 stability multipliers.
%
% PROTOTYPE:
%   result = rpo.compute_stability(cfg, a1_result)
%
% DESCRIPTION:
%   Integrate planar and spatial monodromy matrices for the A1 family and
%   reference RPO. Remove the two multipliers nearest unity before forming
%   spectral radii; classify the reference orbit using the frozen tolerance.
%
% INPUT:
%   cfg         Configuration from reproduction_config(); uses constants
%               and reference_rpo.stability settings.
%   a1_result   A1 result supplying data.family and data.rpo.state / period_tu.
%
% OUTPUT:
%   result      Struct with family multipliers, spectral radii, periods [days],
%               and perilune altitudes [km] in data. Multipliers are [-].
%               metrics.reference_rpo holds radii and stability flags.
%               Configuration and provenance are attached separately.

% Read the family and allocate the nontrivial multiplier arrays.
family = a1_result.data.family;
settings = cfg.reference_rpo.stability;
constants = cfg.constants;
use_parallel = rpo.parallel_available(settings.use_parallel);
n_orbits = size(family, 1);

planar_multipliers = NaN(n_orbits, 2);
spatial_multipliers = NaN(n_orbits, 4);
perilune_altitude_km = NaN(n_orbits, 1);

% Integrate each orbit over its full period to form the monodromy matrices.
fprintf('Linear stability: %d A1 orbits using %s execution.\n', ...
    n_orbits, local_execution_name(use_parallel));
if use_parallel
    parfor index = 1:n_orbits
        [planar_multipliers(index, :), spatial_multipliers(index, :), ...
            perilune_altitude_km(index)] = local_orbit_stability( ...
            family(index, :), constants);
    end
else
    for index = 1:n_orbits
        [planar_multipliers(index, :), spatial_multipliers(index, :), ...
            perilune_altitude_km(index)] = local_orbit_stability( ...
                family(index, :), constants);
    end
end

% Reduce each family spectrum and evaluate the corrected RPO separately.
planar_spectral_radius = max(abs(planar_multipliers), [], 2);
spatial_spectral_radius = max(abs(spatial_multipliers), [], 2);
period_days = 2 * family(:, 3) * cfg.constants.time_unit_days;
reference_state = a1_result.data.rpo.state;
reference_orbit = [reference_state(1), reference_state(4), ...
    a1_result.data.rpo.period_tu / 2, NaN, 0];
[reference_planar, reference_spatial] = local_orbit_stability( ...
    reference_orbit, constants);
reference_planar_radius = max(abs(reference_planar));
reference_spatial_radius = max(abs(reference_spatial));

% Store family indicators and classify the reference RPO with the frozen tolerance.
result.data.family = family;
result.data.perilune_altitude_km = perilune_altitude_km;
result.data.period_days = period_days;
result.data.planar_multipliers = planar_multipliers;
result.data.spatial_multipliers = spatial_multipliers;
result.data.planar_spectral_radius = planar_spectral_radius;
result.data.spatial_spectral_radius = spatial_spectral_radius;
result.data.branch_number = family(:, 5);
result.data.reference_rpo.planar_multipliers = reference_planar;
result.data.reference_rpo.spatial_multipliers = reference_spatial;

result.metrics.reference_rpo.planar_multipliers = ...
    reference_planar;
result.metrics.reference_rpo.spatial_multipliers = ...
    reference_spatial;
result.metrics.reference_rpo.planar_spectral_radius = ...
    reference_planar_radius;
result.metrics.reference_rpo.spatial_spectral_radius = ...
    reference_spatial_radius;
result.metrics.reference_rpo.planar_stable = ...
    reference_planar_radius <= 1 + settings.tolerance;
result.metrics.reference_rpo.spatial_stable = ...
    reference_spatial_radius <= 1 + settings.tolerance;
end

function [planar, spatial, altitude_km] = local_orbit_stability(orbit, constants)
%LOCAL_ORBIT_STABILITY Evaluate nontrivial planar/spatial multipliers.
%
% PROTOTYPE:
%   [planar, spatial, altitude_km] = local_orbit_stability(orbit, constants)
%
% DESCRIPTION:
%   Propagate the same planar orbit in 2D and embedded 3D coordinates.
%   Each full-period STM maps perturbations from the starting epoch back
%   to that epoch after one orbit. Remove the two multipliers nearest +1.
%
% INPUT:
%   orbit      [1 x 5] Row [x0, vy0, half-period, final x, k] in DU, VU, TU;
%                only the first three entries are used.
%   constants  Earth-Moon mass ratio, distance unit [km], and Moon radius [km].
%
% OUTPUT:
%   planar       [1 x 2] Nontrivial planar Floquet multipliers [-].
%   spatial      [1 x 4] Nontrivial spatial Floquet multipliers [-].
%   altitude_km  Half-period x-axis perilune-altitude estimate [km].

mu = constants.mu_earth_moon;

% Evaluate the half-period state and the planar full-period STM.
state_2d = [orbit(1); 0; 0; orbit(2)];
[half_state, ~] = propagate(state_2d, 0, orbit(3), mu, ...
    '2D', 'CR3BP', 'false', false);
[~, planar_stm] = propagate(state_2d, 0, 2 * orbit(3), mu, ...
    '2D', 'CR3BP', 'true', false);
planar = local_nontrivial_multipliers(planar_stm, 2);

% Embed the same orbit in 3D and include its vertical perturbation modes.
state_3d = [orbit(1); 0; 0; 0; orbit(2); 0];
[~, spatial_stm] = propagate(state_3d, 0, 2 * orbit(3), mu, ...
    '3D', 'CR3BP', 'true', false);
spatial = local_nontrivial_multipliers(spatial_stm, 4);

% Retain the frozen x-axis half-period perilune-altitude convention.
altitude_km = (1 - mu - half_state(1)) ...
    * constants.distance_unit_km - constants.moon_radius_km;
end

function multipliers = local_nontrivial_multipliers(monodromy, count)
%LOCAL_NONTRIVIAL_MULTIPLIERS Remove the two eigenvalues nearest +1.
%
% monodromy is the full-period STM. Return count remaining multipliers
% as a row vector; count is 2 for planar and 4 for spatial dynamics.

all_multipliers = eig(monodromy);

% Remove the two autonomous-system unit multipliers by proximity to +1.
[~, order] = sort(abs(all_multipliers - 1), 'ascend');
all_multipliers(order(1:2)) = [];
multipliers = reshape(all_multipliers(1:count), 1, []);
end

function name = local_execution_name(use_parallel)
%LOCAL_EXECUTION_NAME Return the serial or parallel label for progress output.

if use_parallel
    name = 'parallel';
else
    name = 'serial';
end
end
