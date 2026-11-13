function setup = build_vleo_rpo_ms_initial_guess_2d(cfg, seed_package, x0_rpo_planar_3bp)
%BUILD_VLEO_RPO_MS_INITIAL_GUESS_2D Build the planar transfer warm start and mapped NLP setup.
%
% PROTOTYPE:
%   setup = build_vleo_rpo_ms_initial_guess_2d(cfg, seed_package, x0_rpo_planar_3bp)
%
% DESCRIPTION:
%   Shift the patched seed to departure epoch zero, recover the corresponding
%   initial Sun angle, and sample shooting nodes with the selected placement
%   rule. Project the terminal seed onto the free-phase CR3BP RPO and map
%   endpoint velocity jumps into planar Oshima variables. Endpoint states are
%   reconstructed from boundary phases; only interior node states are packed.
%
% INPUT:
%   cfg                   Validation configuration from cfg.earth_to_rpo, with
%                         node, departure, timing, impulse, and solver settings.
%   seed_package          Optional earth2rpo.build_seed output. If omitted,
%                         load seed_package from cfg.seed_path. Rotating state
%                         rows use [x, y, vx, vy] in DU/VU; time vectors use TU.
%   x0_rpo_planar_3bp     Optional [4 elements] canonical CR3BP target state.
%                         If omitted, load xx0_RPO from cfg.rpo_3bp_seed_path.
%
% OUTPUT:
%   setup                 Struct with cfg, seed_package, constants, nlp,
%                         y0/lb/ub, sampled node states/times, phase guesses,
%                         active impulse mask, and node-density diagnostics.
%                         nlp.map defines the packed decision-vector layout;
%                         nlp carries physical scales and constraint policies.
%                         Departure/Sun angles use rad, RPO phase and segment
%                         durations use TU, optional altitude uses km.

if nargin < 2
    seed_data = load(cfg.seed_path);
    seed_package = seed_data.seed_package;
end

% Store the physical and normalization constants.
constants = local_build_constants(seed_package, cfg);

% Shift the saved seed so the departure maneuver starts at t = 0.
seed_time_tu = seed_package.concatenated_time_tu(:);
seed_state_rot = seed_package.concatenated_state_rot;
seed_time_shifted_tu = seed_time_tu - seed_time_tu(1);

% Read the seed time of flight from the shifted timeline.
seed_tof_tu = seed_time_shifted_tu(end);

% Recover the boundary-crossing time from the departure epoch.
match_time_shifted_tu = seed_package.match_time_tu - seed_time_tu(1);

% Rebuild the Sun phase at the departure epoch.
theta_sun_boundary = seed_package.selected_candidate.thetaSun0;
theta_sun0_ig = mod(theta_sun_boundary - constants.ws_sun * match_time_shifted_tu, 2 * pi);

% Load the canonical planar 3BP RPO seed used by the optimizer target.
if nargin < 3
    x0_rpo_planar_3bp = local_load_x0_rpo_planar_3bp(cfg);
end

% Build the seed nodes with the requested placement rule.
[t_nodes_ig, X_nodes_ig, node_density_data] = local_build_seed_nodes(seed_time_shifted_tu, seed_state_rot, cfg, constants);

% Read the inertial periapsis angle saved at the original departure epoch.
phase_departure_inertial_ig = atan2(seed_package.periapsis_state_inertial(2), seed_package.periapsis_state_inertial(1));

% Shift the inertial periapsis angle into the rotating frame used at t = 0.
phase_departure_ig = mod(phase_departure_inertial_ig - seed_time_tu(1), 2 * pi);

% Read the seed departure altitude from the saved selected candidate.
departure_altitude_seed_km = seed_package.selected_candidate.periapsis_altitude_km;

% Select the departure altitude used by the active warm start.
departure_altitude_ig_km = local_select_departure_altitude_ig(cfg, departure_altitude_seed_km);

% Store the requested fixed target altitude for modes that enforce one.
departure_altitude_target_km = cfg.departure_altitude_km;

% Project the saved terminal state onto the free-phase 3BP RPO.
[phase_rpo_ig, x_rpo_target_ig] = local_project_seed_terminal_to_rpo_3bp(X_nodes_ig(end, :).', x0_rpo_planar_3bp, constants, cfg);

% Build the Earth departure state that corresponds to the warm start.
[x_departure_pre_ig, ~, ~] = circular_parking_state_2d(phase_departure_ig, departure_altitude_ig_km, constants);

% Convert the seed departure velocity jump into planar Oshima variables.
dv_departure_ig = X_nodes_ig(1, 3:4).' - x_departure_pre_ig(3:4);
[U_departure_ig, is_case2b_departure] = dv_to_planar_oshima_u_2d(dv_departure_ig);

% Convert the seed arrival insertion jump into planar Oshima variables.
dv_arrival_ig = x_rpo_target_ig(3:4) - X_nodes_ig(end, 3:4).';
[U_arrival_ig, is_case2b_arrival] = dv_to_planar_oshima_u_2d(dv_arrival_ig);

% Initialize all node-wise regularized variables near zero.
U_nodes_ig = eps * ones(cfg.n_nodes, 2);
U_nodes_ig(1, :) = U_departure_ig(:).';
U_nodes_ig(end, :) = U_arrival_ig(:).';

% Keep the initial case-switch set fixed during the solve.
I_case2b = zeros(2, 1);
case2b_count = 0;
if is_case2b_departure
    case2b_count = case2b_count + 1;
    I_case2b(case2b_count) = 1;
end
if is_case2b_arrival
    case2b_count = case2b_count + 1;
    I_case2b(case2b_count) = cfg.n_nodes;
end
I_case2b = I_case2b(1:case2b_count);

% Activate all node impulses or only the endpoint impulses.
active_impulse_mask = false(cfg.n_nodes, 1);
if strcmpi(cfg.impulse_mode, 'every_node')
    active_impulse_mask(:) = true;
else
    active_impulse_mask([1, cfg.n_nodes]) = true;
end

% Pack the warm start into the NLP decision vector.
[y0, map] = local_pack_decision_vector(X_nodes_ig, U_nodes_ig, t_nodes_ig, phase_departure_ig, phase_rpo_ig, theta_sun0_ig, departure_altitude_ig_km, cfg);

% Build the lower and upper bounds for the NLP variables.
[lb, ub] = local_build_bounds(y0, map, cfg, constants, departure_altitude_seed_km, departure_altitude_target_km);

% Store the full NLP description used by the optimizer.
nlp = struct;
nlp.N = cfg.n_nodes;
nlp.map = map;
nlp.mu = constants.mu;
nlp.muE = constants.muE;
nlp.DU = constants.DU;
nlp.VU = constants.VU;
nlp.R_e = constants.R_e;
nlp.R_m = constants.R_m;
nlp.T_RPO = constants.T_RPO;
nlp.T_RPO_3bp = constants.T_RPO_3bp;
nlp.ws_sun = constants.ws_sun;
nlp.abs_tol = cfg.abs_tol;
nlp.rel_tol = cfg.rel_tol;
nlp.min_dt_tu = cfg.min_dt_tu;
nlp.seed_tof_tu = seed_tof_tu;
nlp.max_tof_tu = cfg.tof_max_scale * seed_tof_tu;
nlp.moon_min_altitude_km = cfg.moon_min_altitude_km;
nlp.rmin_moon_du = (constants.R_m + cfg.moon_min_altitude_km) / constants.DU;
nlp.departure_constraint_mode = cfg.departure_constraint_mode;
nlp.departure_altitude_km_fixed = local_select_departure_altitude_fixed(cfg, departure_altitude_seed_km);
nlp.departure_altitude_target_km = departure_altitude_target_km;
nlp.active_impulse_mask = active_impulse_mask;
nlp.I_case2b = I_case2b(:);
nlp.impulse_mode = cfg.impulse_mode;
nlp.tof_mode = cfg.tof_mode;
nlp.x0_rpo_planar_3bp = x0_rpo_planar_3bp(:);

% Package the complete setup used by the root script.
setup = struct;
setup.cfg = cfg;
setup.seed_package = seed_package;
setup.constants = constants;
setup.nlp = nlp;
setup.y0 = y0;
setup.lb = lb;
setup.ub = ub;
setup.seed_time_shifted_tu = seed_time_shifted_tu;
setup.seed_state_rot = seed_state_rot;
setup.seed_match_time_shifted_tu = match_time_shifted_tu;
setup.seed_tof_tu = seed_tof_tu;
setup.X_nodes_ig = X_nodes_ig;
setup.t_nodes_ig = t_nodes_ig;
setup.node_density_data = node_density_data;
setup.theta_sun0_ig = theta_sun0_ig;
setup.phase_departure_ig = phase_departure_ig;
setup.phase_rpo_ig = phase_rpo_ig;
setup.departure_altitude_seed_km = departure_altitude_seed_km;
setup.departure_altitude_ig_km = departure_altitude_ig_km;
setup.departure_altitude_target_km = departure_altitude_target_km;
setup.active_impulse_mask = active_impulse_mask;
setup.x0_rpo_planar_3bp = x0_rpo_planar_3bp(:);

end

function constants = local_build_constants(seed_package, cfg)
%LOCAL_BUILD_CONSTANTS Collect the seed's physical scales and the configured CR3BP target period.
%
% DESCRIPTION:
%   Preserve DU [km], VU [km/s], TU_days [days/TU], body radii [km],
% normalized mass ratio, solar angular rate [rad/TU], and RPO periods [TU].

constants = struct;
constants.mu = seed_package.constants.mu;
constants.DU = seed_package.constants.DU;
constants.VU = seed_package.constants.VU;
constants.TU_days = seed_package.constants.TU_days;
constants.R_e = seed_package.constants.R_e;
constants.R_m = seed_package.constants.R_m;
constants.T_RPO = seed_package.constants.T_RPO;
constants.T_RPO_3bp = cfg.T_RPO_3bp_tu;

% Store the Earth gravitational parameter used by the parking orbit.
constants.muE = 398600;

% Store the nondimensional Sun angular velocity used by the BCR4BP model.
constants.ws_sun = -9.25195985e-1;

end

function x0_rpo_planar_3bp = local_load_x0_rpo_planar_3bp(cfg)
%LOCAL_LOAD_X0_RPO_PLANAR_3BP Load the canonical planar target from the configured MAT file.
%
% DESCRIPTION:
%   Read xx0_RPO as [x; y; vx; vy] in rotating barycentric DU/VU.

seed_data = load(cfg.rpo_3bp_seed_path);

% Read the planar seed from the expected field.
x0_rpo_planar_3bp = seed_data.xx0_RPO(:);

end

function [phase_rpo_ig, x_rpo_target_ig] = local_project_seed_terminal_to_rpo_3bp(x_seed_terminal, x0_rpo_planar_3bp, constants, cfg)
%LOCAL_PROJECT_SEED_TERMINAL_TO_RPO_3BP Choose a CR3BP target phase near the terminal patched seed state.
%
% PROTOTYPE:
%   [phase_rpo_ig, x_rpo_target_ig] = local_project_seed_terminal_to_rpo_3bp( ...
%       x_seed_terminal, x0_rpo_planar_3bp, constants, cfg)
%
% DESCRIPTION:
%   Sample one nominal target period, minimize the unweighted full-state
%   mismatch on that sampled orbit, wrap the selected phase, then repropagate
%   exactly to that selected epoch. This is a sampled warm-start projection,
%   not an additional phase optimization.
%
% INPUT:
%   x_seed_terminal       [4 elements] Terminal rotating seed [DU/VU].
%   x0_rpo_planar_3bp     [4 elements] Canonical CR3BP RPO initial state [DU/VU].
%   constants             Mass ratio and T_RPO_3bp [TU].
%   cfg                   abs_tol/rel_tol integration settings.
%
% OUTPUT:
%   phase_rpo_ig          Selected target time offset within one period [TU].
%   x_rpo_target_ig       [4 x 1] Target state at that phase [DU/VU].

[~, ~, ~, x_rpo_guess, t_rpo_guess] = propagate(x0_rpo_planar_3bp, 0, constants.T_RPO_3bp, constants.mu, '2D', 'CR3BP', 'false', false, [], [], cfg.abs_tol, cfg.rel_tol);

% Measure the full-state mismatch between the seed end state and the 3BP orbit.
state_error = x_rpo_guess - x_seed_terminal(:).';
state_error_norm = sqrt(sum(state_error .^ 2, 2));

% Read the nearest 3BP phase sample from the propagated orbit.
[~, best_idx] = min(state_error_norm);
phase_rpo_ig = t_rpo_guess(best_idx);
phase_rpo_ig = mod(phase_rpo_ig, constants.T_RPO_3bp);

% Re-evaluate the 3BP target state exactly at the selected free phase.
x_rpo_target_ig = propagate(x0_rpo_planar_3bp, 0, phase_rpo_ig, constants.mu, '2D', 'CR3BP', 'false', false, [], [], cfg.abs_tol, cfg.rel_tol);

end

function [t_nodes_ig, X_nodes_ig, node_density_data] = local_build_seed_nodes(seed_time_shifted_tu, seed_state_rot, cfg, constants)
%LOCAL_BUILD_SEED_NODES Sample the seed using uniform-time or body-proximity node placement.
%
% PROTOTYPE:
%   [t_nodes_ig, X_nodes_ig, node_density_data] = local_build_seed_nodes( ...
%       seed_time_shifted_tu, seed_state_rot, cfg, constants)
%
% DESCRIPTION:
%   For body-proximity placement, integrate the positive Earth/Moon proximity
%   density and sample uniformly in its cumulative coordinate. Interpolate
%   the original seed states at the resulting epochs with pchip.
%
% INPUT:
%   seed_time_shifted_tu  [Nt x 1] Seed timeline [TU], beginning at zero.
%   seed_state_rot        [Nt x 4] Matching rotating states [x, y, vx, vy].
%   cfg                   Node count, placement mode, density scales/weights.
%   constants             mu and DU [km], for physical primary distances.
%
% OUTPUT:
%   t_nodes_ig           [N x 1] Sampled node epochs [TU], retaining endpoints.
%   X_nodes_ig           [N x 4] Interpolated rotating states in DU/VU.
%   node_density_data    Seed-time/proximity/density/cumulative-coordinate
%                          arrays used for placement diagnostics.

seed_tof_tu = seed_time_shifted_tu(end);

% Build the uniform-time nodes when requested.
if strcmpi(cfg.node_placement, 'uniform_time')
    % Place all nodes uniformly in the shifted seed time.
    t_nodes_ig = linspace(0, seed_tof_tu, cfg.n_nodes).';

    % Store a flat density profile for reporting and plots.
    earth_proximity = zeros(size(seed_time_shifted_tu));
    moon_proximity = zeros(size(seed_time_shifted_tu));
    density_profile = ones(size(seed_time_shifted_tu));
    density_coordinate = seed_time_shifted_tu / seed_tof_tu;

else
    % Read the Earth-centered seed distance in kilometers.
    earth_distance_km = constants.DU * hypot(seed_state_rot(:, 1) + constants.mu, seed_state_rot(:, 2));

    % Read the Moon-centered seed distance in kilometers.
    moon_distance_km = constants.DU * hypot(seed_state_rot(:, 1) - 1 + constants.mu, seed_state_rot(:, 2));

    % Build the smooth Earth-proximity activation used for node clustering.
    earth_proximity = 1 ./ (1 + (earth_distance_km / cfg.node_density_earth_scale_km) .^ cfg.node_density_power);

    % Build the smooth Moon-proximity activation used for node clustering.
    moon_proximity = 1 ./ (1 + (moon_distance_km / cfg.node_density_moon_scale_km) .^ cfg.node_density_power);

    % Build the positive density profile used to remap the node times.
    density_profile = 1 + cfg.node_density_strength * (cfg.node_density_earth_weight * earth_proximity + cfg.node_density_moon_weight * moon_proximity);

    % Integrate the density profile into a monotone placement coordinate.
    density_coordinate = cumtrapz(seed_time_shifted_tu, density_profile);
    density_coordinate = density_coordinate / density_coordinate(end);

    % Place the nodes uniformly in the density coordinate.
    t_nodes_ig = interp1(density_coordinate, seed_time_shifted_tu, linspace(0, 1, cfg.n_nodes).', 'linear');
    t_nodes_ig(1) = 0;
    t_nodes_ig(end) = seed_tof_tu;
end

% Sample the seed state on the requested node times.
X_nodes_ig = interp1(seed_time_shifted_tu, seed_state_rot, t_nodes_ig, 'pchip');

% Store the node-density diagnostics in one compact structure.
node_density_data = struct;
node_density_data.time_tu = seed_time_shifted_tu;
node_density_data.earth_proximity = earth_proximity;
node_density_data.moon_proximity = moon_proximity;
node_density_data.density_profile = density_profile;
node_density_data.density_coordinate = density_coordinate;

end

function departure_altitude_km = local_select_departure_altitude_ig(cfg, departure_altitude_seed_km)
%LOCAL_SELECT_DEPARTURE_ALTITUDE_IG Choose the departure-altitude warm start from the active mode.
%
% DESCRIPTION:
%   Use the seed, its clipped interval value, or configured fixed altitude [km].

if strcmpi(cfg.departure_constraint_mode, 'seed_altitude') || strcmpi(cfg.departure_constraint_mode, 'fixed_altitude_seed_warmstart')
    departure_altitude_km = departure_altitude_seed_km;

% Clip the seed altitude to the permitted interval when altitude is free.
elseif strcmpi(cfg.departure_constraint_mode, 'altitude_interval')
    departure_altitude_km = min(max(departure_altitude_seed_km, cfg.departure_altitude_bounds_km(1)), cfg.departure_altitude_bounds_km(2));

% Use the configured fixed altitude otherwise.
else
    departure_altitude_km = cfg.departure_altitude_km;
end

end

function departure_altitude_km_fixed = local_select_departure_altitude_fixed(cfg, departure_altitude_seed_km)
%LOCAL_SELECT_DEPARTURE_ALTITUDE_FIXED Choose the fixed departure altitude [km] used without an altitude variable.

if strcmpi(cfg.departure_constraint_mode, 'seed_altitude')
    departure_altitude_km_fixed = departure_altitude_seed_km;

% Use the configured fixed altitude otherwise.
else
    departure_altitude_km_fixed = cfg.departure_altitude_km;
end

end

function [y0, map] = local_pack_decision_vector(X_nodes_ig, U_nodes_ig, t_nodes_ig, phase_departure_ig, phase_rpo_ig, theta_sun0_ig, departure_altitude_ig_km, cfg)
%LOCAL_PACK_DECISION_VECTOR Pack endpoint impulses, interior states, durations, and boundary phases.
%
% PROTOTYPE:
%   [y0, map] = local_pack_decision_vector(X_nodes_ig, U_nodes_ig, t_nodes_ig, ...
%       phase_departure_ig, phase_rpo_ig, theta_sun0_ig, ...
%       departure_altitude_ig_km, cfg)
%
% DESCRIPTION:
%   Record each variable index while building the warm start. Endpoint states
%   are not independent variables: departure comes from its circular phase,
%   arrival from CR3BP RPO phase. Include departure altitude only for the
%   two supported variable-altitude modes.
%
% INPUT:
%   X_nodes_ig              [N x 4] Rotating seed states [DU/VU].
%   U_nodes_ig              [N x 2] Regularized impulse coefficients.
%   t_nodes_ig              [N x 1] Seed epochs [TU].
%   phase_departure_ig      Departure-circle angle [rad].
%   phase_rpo_ig            CR3BP target time offset [TU].
%   theta_sun0_ig           Initial Sun angle [rad].
%   departure_altitude_ig_km Departure warm-start altitude [km].
%   cfg                     Departure constraint mode.
%
% OUTPUT:
%   y0    Packed initial decision-vector column.
%   map   Struct with node(j).X/U, dt_nodes, phase_departure, phase_rpo,
%           theta_sun0, and optional departure_altitude_km indices.
%           Endpoint node(j).X arrays are empty.

N = size(X_nodes_ig, 1);

% Preallocate the map structure.
map = struct;
map.node = repmat(struct('X', [], 'U', []), N, 1);
offset = 0;

% Pack the first node regularized variables.
idx_u1 = offset + (1:2);
map.node(1).U = idx_u1;
offset = idx_u1(end);

% Pack each interior state and its regularized variables.
for j = 2:(N - 1)
    % Store the interior node state.
    idx_x = offset + (1:4);
    map.node(j).X = idx_x;
    offset = idx_x(end);

    % Store the interior node regularized variables.
    idx_u = offset + (1:2);
    map.node(j).U = idx_u;
    offset = idx_u(end);
end

% Pack the last node regularized variables.
idx_un = offset + (1:2);
map.node(N).U = idx_un;
offset = idx_un(end);

% Pack the segment durations.
idx_dt = offset + (1:(N - 1));
map.dt_nodes = idx_dt;
offset = idx_dt(end);

% Pack the scalar phase variables.
map.phase_departure = offset + 1;
map.phase_rpo = offset + 2;
map.theta_sun0 = offset + 3;
offset = map.theta_sun0;

% Pack the departure altitude only when it is a decision variable.
if strcmpi(cfg.departure_constraint_mode, 'altitude_interval') || strcmpi(cfg.departure_constraint_mode, 'fixed_altitude_seed_warmstart')
    map.departure_altitude_km = offset + 1;
    offset = map.departure_altitude_km;
else
    map.departure_altitude_km = [];
end

% Allocate the full decision vector.
y0 = zeros(offset, 1);

% Store the first node regularized variables.
y0(map.node(1).U) = U_nodes_ig(1, :).';

% Store the interior node states and regularized variables.
for j = 2:(N - 1)
    % Store the interior state at node j.
    y0(map.node(j).X) = X_nodes_ig(j, :).';

    % Store the regularized variables at node j.
    y0(map.node(j).U) = U_nodes_ig(j, :).';
end

% Store the last node regularized variables.
y0(map.node(N).U) = U_nodes_ig(N, :).';

% Store the segment durations.
y0(map.dt_nodes) = diff(t_nodes_ig(:));

% Store the scalar phase variables.
y0(map.phase_departure) = phase_departure_ig;
y0(map.phase_rpo) = phase_rpo_ig;
y0(map.theta_sun0) = theta_sun0_ig;

% Store the departure altitude when it is a decision variable.
if ~isempty(map.departure_altitude_km)
    y0(map.departure_altitude_km) = departure_altitude_ig_km;
end

end

function [lb, ub] = local_build_bounds(y0, map, cfg, constants, departure_altitude_seed_km, departure_altitude_target_km)
%LOCAL_BUILD_BOUNDS Build duration, boundary-phase, and optional altitude box bounds.
%
% PROTOTYPE:
%   [lb, ub] = local_build_bounds(y0, map, cfg, constants, ...
%       departure_altitude_seed_km, departure_altitude_target_km)
%
% DESCRIPTION:
%   Keep segment durations positive, departure/Sun angles within one turn,
%   and target phase within one nominal period. For altitude warm starts,
%   include both seed/target altitudes with the configured margin; otherwise
%   use the explicitly configured free-altitude interval.
%
% INPUT:
%   y0 / map        Packed seed and its variable-index map.
%   cfg             Minimum duration and departure-altitude policy/settings.
%   constants       Nominal CR3BP RPO period [TU].
%   departure_altitude_seed_km / departure_altitude_target_km
%                   Seed/target altitudes [km] used by the warm-start mode.
%
% OUTPUT:
%   lb, ub          [numel(y0) x 1] Bounds in each variable's native units.
%                     State and impulse variables remain otherwise unbounded.

lb = -inf(size(y0));
ub = inf(size(y0));

% Keep the segment durations strictly positive.
lb(map.dt_nodes) = cfg.min_dt_tu;

% Keep the departure phase inside one full revolution.
lb(map.phase_departure) = 0;
ub(map.phase_departure) = 2 * pi;

% Keep the RPO phase inside one reference period.
lb(map.phase_rpo) = 0;
ub(map.phase_rpo) = constants.T_RPO_3bp;

% Keep the initial Sun phase inside one full revolution.
lb(map.theta_sun0) = 0;
ub(map.theta_sun0) = 2 * pi;

% Bound the free departure altitude when it is active.
if ~isempty(map.departure_altitude_km)
    if strcmpi(cfg.departure_constraint_mode, 'fixed_altitude_seed_warmstart')
        % Build automatic bounds that contain both the seed and target altitudes.
        altitude_min_km = min(departure_altitude_seed_km, departure_altitude_target_km);
        altitude_max_km = max(departure_altitude_seed_km, departure_altitude_target_km);
        lb(map.departure_altitude_km) = max(0, altitude_min_km - cfg.departure_altitude_warmstart_margin_km);
        ub(map.departure_altitude_km) = altitude_max_km + cfg.departure_altitude_warmstart_margin_km;
    else
        % Use the user-provided interval for genuinely free-altitude solves.
        lb(map.departure_altitude_km) = cfg.departure_altitude_bounds_km(1);
        ub(map.departure_altitude_km) = cfg.departure_altitude_bounds_km(2);
    end
end

end
