function problem_definition = build_casadi_rpo2ds_initial_state_problem(y0, lb, ub, problem, opts)
%BUILD_CASADI_RPO2DS_INITIAL_STATE_PROBLEM Build the symbolic CR3BP multiple-shooting NLP.
%
% PROTOTYPE:
%   problem_definition = build_casadi_rpo2ds_initial_state_problem(y0, lb, ub, ...
%       problem, opts)
%
% DESCRIPTION:
%   Decode the existing node/duration layout, build normalized CR3BP flow
%   maps, and preserve the shared constraint ordering. CasADi differentiates
%   the flows and algebraic boundary constraints. The regularized opt cost
%   sums impulse magnitudes; reach mode maximizes final v-infinity.
%
% INPUT:
%   y0, lb, ub  [nvars x 1] Seed and variable bounds from the transfer setup.
%   problem     Transfer problem_data struct; first 9*n entries are [X(6); U(3)],
%               dt_idx selects durations [TU], and time/phase/target fields
%               define the departure and arrival boundaries. See
%               RPO2DS_transfer and constr_CR3BP_phase_core for conventions.
%   opts        Solver options; optional casadi_path selects the interface.
%
% OUTPUT:
%   problem_definition  Struct with nlp, lbx/ubx, lbg/ubg, objective and
%                         constraint evaluator functions, and dimensions.
%                         Constraint bounds encode ceq = 0 and c <= 0;
%                         this function builds the NLP without solving it.

if isfield(opts, 'casadi_path') && ~isempty(opts.casadi_path)
    addpath(opts.casadi_path);
end
import casadi.*

% Read the fixed NLP dimensions from the existing RPO2DS setup package.
n = problem.n;
nvars = numel(y0);
dt_idx = problem.dt_idx(:);
t1_idx = problem.t1_idx;
time_anchor_mode = 'initial';
if isfield(problem, 'time_anchor_mode') && ~isempty(problem.time_anchor_mode)
    time_anchor_mode = char(lower(string(problem.time_anchor_mode)));
end
has_final_time_lock = local_has_final_time_lock(problem, time_anchor_mode);

% Create one symbolic decision vector with the existing packed layout.
y = MX.sym('y', nvars, 1);
node_blocks = reshape(y(1:9 * n), 9, n);
X_nodes = node_blocks(1:6, :);
U_nodes = node_blocks(7:9, :);
dt_nodes = y(dt_idx);
if strcmp(time_anchor_mode, 'final')
    t1 = problem.final_epoch_time - sum1(dt_nodes);
    t_nodes = t1 + [0; cumsum(dt_nodes)];
else
    t1 = y(t1_idx);
    t_nodes = t1 + [0; cumsum(dt_nodes)];
end
active_impulse_mask = local_active_impulse_mask(problem, n);

% Build the two CasADi flow maps used by departure and segment constraints.
integrator_options = local_build_integrator_options();
rpo_flow = local_build_cr3bp_flow('rpo2ds_rpo_flow', problem, integrator_options);
segment_flow = local_build_cr3bp_flow('rpo2ds_segment_flow', problem, integrator_options);

% Build the regularized delta-v objective used by the final optimization stage.
objective = 0;
for j = 1:n
    if active_impulse_mask(j)
        objective = objective + sumsqr(U_nodes(:, j));
    end
end

% Switch the final-stage objective to maximum reachable v-infinity in reach mode.
if strcmp(problem.mode, 'reach')
    objective = -y(problem.vinf_mag_idx);
end

% Build the departure state from the active fixed, free, or time-derived phase source.
x_departure = local_departure_state_symbolic(y, t1, problem, rpo_flow);

% Store position defects separately so their order matches fmincon exactly.
segment_position_defects = cell(n - 1, 1);
segment_velocity_defects = cell(n - 1, 1);
for i = 1:n - 1
    % Propagate one multiple-shooting segment over its optimized duration.
    segment_result = segment_flow('x0', X_nodes(:, i), 'p', dt_nodes(i));
    x_propagated = segment_result.xf;

    % Apply the regularized impulse at the downstream node.
    dv_next = local_regularized_impulse_if_active_symbolic(U_nodes(:, i + 1), ...
        any(problem.I_2b == (i + 1)), active_impulse_mask(i + 1));
    segment_position_defects{i} = x_propagated(1:3) - X_nodes(1:3, i + 1);
    segment_velocity_defects{i} = X_nodes(4:6, i + 1) - x_propagated(4:6) - dv_next;
end

% Match the first node to the phase-derived RPO departure state plus burn.
dv_departure = local_regularized_impulse_if_active_symbolic(U_nodes(:, 1), any(problem.I_2b == 1), active_impulse_mask(1));
departure_defect = [X_nodes(1:3, 1) - x_departure(1:3); ...
    X_nodes(4:6, 1) - x_departure(4:6) - dv_departure];

% Build the Earth-periapsis and outgoing v-infinity constraints.
Xn = X_nodes(:, n);
tn = t_nodes(end);
[r1, r1dot, r1ddot] = local_primary_radius_derivatives_symbolic(Xn, problem.mu, -problem.mu);
has_final_earth_distance_cap = local_has_standard_fixed_phase_time_final_earth_cap(problem);
Xn_in = local_rot2in_symbolic(Xn, tn, problem.mu);
rp_earth = Xn_in(1:3) * problem.DU;
vp_earth = Xn_in(4:6) * problem.VU;
vinf_out = local_earth_vinf_out_symbolic(rp_earth, vp_earth, problem.muE);
vinf_target = local_active_vinf_target_symbolic(y, t_nodes(end), problem);
vinf_defect = vinf_out - vinf_target;

% Build the Moon-periapsis constraints on the selected fixed node.
X_moon = X_nodes(:, problem.idx_peri_moon);
[r2, r2dot, r2ddot] = local_primary_radius_derivatives_symbolic(X_moon, problem.mu, 1 - problem.mu);

% Keep the equality vector order identical to constr_CR3BP_phase_core.
ceq = [vertcat(segment_position_defects{:}); ...
    vertcat(segment_velocity_defects{:}); ...
    departure_defect; ...
    r1dot; ...
    vinf_defect; ...
    r2dot];

% Add the specified final epoch row after the shared event constraints.
if has_final_time_lock
    ceq = [ceq; t_nodes(end) - problem.t_nodeN_fixed];
end

% Keep the base inequality vector order identical to constr_CR3BP_phase_core.
c_base = [-dt_nodes; ...
    -r1ddot; ...
    problem.rminEarth^2 - r1^2; ...
    sum(dt_nodes) - problem.maxTransferTOF; ...
    problem.rminMoon^2 - r2^2; ...
    -r2ddot];

% Start the inequality vector with the shared base rows.
c = c_base;

% Add the IPOPT-only terminal Earth-distance cap for the fixed phase/time profile.
if has_final_earth_distance_cap
    c = [c; r1^2 - problem.maxFinalEarthDistance^2];
end

% Add the operational initial-state timing rows only for that profile.
if local_has_initial_state_timing(problem)
    c = [c; ...
        t_nodes(end) - problem.initial_epoch_time - problem.maxTotalTime; ...
        problem.initial_epoch_time - t1];
end

% Add the RPO wait-only cap when the initial-state RTH profile requests it.
if local_has_initial_wait_cap(problem)
    c = [c; t1 - problem.initial_epoch_time - problem.initial_wait_cap];
end

% Add the regularized maneuver-budget constraint used by reachability mode.
if strcmp(problem.mode, 'reach')
    regularized_budget = 0;
    for j = 1:n
        regularized_budget = regularized_budget + sumsqr(U_nodes(:, j));
    end
    c = [c; regularized_budget - problem.dv_budget];
end

% Package the full CasADi NLP and helper evaluators for tests.
constraints = [ceq; c];
n_equalities = 6 * n + 5 + has_final_time_lock;
n_inequalities = (n + 4) + has_final_earth_distance_cap + ...
    2 * local_has_initial_state_timing(problem) + local_has_initial_wait_cap(problem) + ...
    strcmp(problem.mode, 'reach');
nlp = struct('x', y, 'f', objective, 'g', constraints);
problem_definition = struct();
problem_definition.nlp = nlp;
problem_definition.lbx = lb(:);
problem_definition.ubx = ub(:);
problem_definition.lbg = [zeros(n_equalities, 1); -inf(n_inequalities, 1)];
problem_definition.ubg = zeros(n_equalities + n_inequalities, 1);
problem_definition.num_equalities = n_equalities;
problem_definition.num_inequalities = n_inequalities;
problem_definition.objective_function = Function('rpo2ds_casadi_objective', {y}, {objective});
problem_definition.constraint_function = Function('rpo2ds_casadi_constraints', {y}, {constraints});
problem_definition.constraint_split_function = Function('rpo2ds_casadi_constraint_split', {y}, {c, ceq});
problem_definition.solver_configuration = struct('integrator_name', 'cvodes', ...
    'integrator_options', integrator_options, ...
    'nvars', nvars, ...
    'n_equalities', n_equalities, ...
    'n_inequalities', n_inequalities);
end

function integrator_options = local_build_integrator_options()
%LOCAL_BUILD_INTEGRATOR_OPTIONS Build the fixed CVODES options used by the IPOPT v1 profile.

integrator_options = struct();
integrator_options.abstol = 1e-12;
integrator_options.reltol = 1e-12;
integrator_options.max_num_steps = 100000;
integrator_options.expand = false;
integrator_options.error_on_fail = true;
integrator_options.print_stats = false;
integrator_options.verbose = false;
end

function x_departure = local_departure_state_symbolic(y, t1, problem, rpo_flow)
%LOCAL_DEPARTURE_STATE_SYMBOLIC Build the symbolic departure state for the active IPOPT profile.

% The fixed-phase standard profile reuses the already phased seed state.
if ~isfield(problem, 'phase_source') || strcmp(problem.phase_source, 'fixed')
    x_departure = problem.x_ini_fixed(:);
    return
end

% The operational initial-state profile derives RPO phase from the first node time.
if strcmp(problem.phase_source, 'initial_time')
    theta_raw = problem.initial_rpo_phase + t1 - problem.initial_epoch_time;
    theta_departure = local_wrap_phase_symbolic(theta_raw, problem);
else
    % The explicit-phase profiles store theta_rpo as one bounded decision variable.
    theta_departure = local_wrap_phase_symbolic(y(problem.phase_theta_idx), problem);
end

% Propagate the planar RPO seed to the selected symbolic departure phase.
departure_result = rpo_flow('x0', problem.x0_rpo_base(:), 'p', theta_departure);
x_departure = departure_result.xf;
end

function theta_wrapped = local_wrap_phase_symbolic(theta_raw, problem)
%LOCAL_WRAP_PHASE_SYMBOLIC Wrap one symbolic RPO phase onto the database geometry interval.

% Keep fixed-phase profiles unchanged when no period metadata is available.
theta_wrapped = theta_raw;
if ~isfield(problem, 'phase_period') || isempty(problem.phase_period) || ...
        ~isfinite(problem.phase_period) || problem.phase_period <= 0
    return
end

% Remove whole-period waits with a differentiable slope away from boundaries.
theta_shifted = theta_raw - problem.phase_origin;
theta_wrapped = problem.phase_origin + theta_shifted - problem.phase_period * floor(theta_shifted / problem.phase_period);
end

function has_timing = local_has_initial_state_timing(problem)
%LOCAL_HAS_INITIAL_STATE_TIMING Check whether the operational initial-state rows are active.

has_timing = isfield(problem, 'phase_source') && ...
    strcmp(problem.phase_source, 'initial_time') && strcmp(problem.operational_mode, 'initial_state');
end

function has_cap = local_has_initial_wait_cap(problem)
%LOCAL_HAS_INITIAL_WAIT_CAP Check whether the initial-state wait-only cap is active.

has_cap = isfield(problem, 'phase_source') && strcmp(problem.phase_source, 'initial_time') && ...
    strcmp(problem.operational_mode, 'initial_state') && isfield(problem, 'initial_wait_cap') && ...
    ~isempty(problem.initial_wait_cap);
end

function has_lock = local_has_final_time_lock(problem, time_anchor_mode)
%LOCAL_HAS_FINAL_TIME_LOCK Check whether propagation times need one explicit final-epoch equality.

has_lock = isfield(problem, 'time_variable_mode') && strcmp(problem.time_variable_mode, 'propagation_times') && ...
    isfield(problem, 'lock_t_nodeN') && logical(problem.lock_t_nodeN) && ...
    ~strcmp(time_anchor_mode, 'final') && isfield(problem, 't_nodeN_fixed') && ~isempty(problem.t_nodeN_fixed);
end

function has_cap = local_has_standard_fixed_phase_time_final_earth_cap(problem)
%LOCAL_HAS_STANDARD_FIXED_PHASE_TIME_FINAL_EARTH_CAP Check whether the IPOPT-only terminal Earth-distance cap should be active.

has_cap = isfield(problem, 'maxFinalEarthDistance') && ~isempty(problem.maxFinalEarthDistance) && ...
    isfinite(problem.maxFinalEarthDistance) && problem.maxFinalEarthDistance > 0 && ...
    isfield(problem, 'phase_source') && strcmp(problem.operational_mode, 'standard') && strcmp(problem.phase_source, 'variable') && ...
    strcmp(problem.geometry_mode, 'prescribed_geometry');
end

function active_impulse_mask = local_active_impulse_mask(problem, n)
%LOCAL_ACTIVE_IMPULSE_MASK Read the active impulse mask and default to every node.

active_impulse_mask = true(n, 1);
if isfield(problem, 'active_impulse_mask') && ~isempty(problem.active_impulse_mask)
    active_impulse_mask = logical(problem.active_impulse_mask(:));
end
end

function flow = local_build_cr3bp_flow(name, problem, integrator_options)
%LOCAL_BUILD_CR3BP_FLOW Build one normalized-duration CR3BP flow map.

import casadi.*

% Scale the autonomous CR3BP dynamics by the symbolic duration parameter.
x = MX.sym([name '_x'], 6, 1);
duration = MX.sym([name '_duration'], 1, 1);
dae = struct();
dae.x = x;
dae.p = duration;
dae.ode = duration * local_rhs_cr3bp_symbolic(x, problem.mu);
flow = integrator(name, 'cvodes', dae, 0, 1, integrator_options);
end

function dxdt = local_rhs_cr3bp_symbolic(X, mu)
%LOCAL_RHS_CR3BP_SYMBOLIC Build the symbolic 3D CR3BP right-hand side.

x = X(1);
y = X(2);
z = X(3);
vx = X(4);
vy = X(5);
vz = X(6);
r1 = sqrt((x + mu)^2 + y^2 + z^2);
r2 = sqrt((x + mu - 1)^2 + y^2 + z^2);
dUdx = x - (1 - mu) * (mu + x) / r1^3 + mu * (1 - mu - x) / r2^3;
dUdy = y - (1 - mu) * y / r1^3 - mu * y / r2^3;
dUdz = -(1 - mu) * z / r1^3 - mu * z / r2^3;
dxdt = [vx; vy; vz; dUdx + 2 * vy; dUdy - 2 * vx; dUdz];
end

function dv_vec = local_regularized_impulse_symbolic(U, is_case2b)
%LOCAL_REGULARIZED_IMPULSE_SYMBOLIC Map one Oshima regularized impulse vector into physical delta-v.

u = U(1);
w = U(2);
s = U(3);
if ~is_case2b
    dv_vec = [u^2 - w^2 - s^2; 2 * u * w; 2 * u * s];
else
    dv_vec = [2 * u * w; u^2 - w^2 - s^2; 2 * u * s];
end
end

function dv_vec = local_regularized_impulse_if_active_symbolic(U, is_case2b, is_active)
%LOCAL_REGULARIZED_IMPULSE_IF_ACTIVE_SYMBOLIC Map an active impulse and return a symbolic zero vector for inactive nodes.

if is_active
    dv_vec = local_regularized_impulse_symbolic(U, is_case2b);
else
    dv_vec = 0 * U;
end
end

function [r, rdot, rddot] = local_primary_radius_derivatives_symbolic(X, mu, primary_x)
%LOCAL_PRIMARY_RADIUS_DERIVATIVES_SYMBOLIC Build symbolic distance, radial velocity, and radial acceleration.

dr = [X(1) - primary_x; X(2); X(3)];
vel = X(4:6);
acc = local_rhs_cr3bp_symbolic(X, mu);
acc = acc(4:6);
r = sqrt(sumsqr(dr));
rv = dr.' * vel;
rdot = rv / r;
rddot = (vel.' * vel + dr.' * acc) / r - rv * rdot / r^2;
end

function X_in = local_rot2in_symbolic(X, t, mu)
%LOCAL_ROT2IN_SYMBOLIC Transform one rotating-frame state into the Earth-centered inertial frame.

c = cos(t);
s = sin(t);
x = X(1);
y = X(2);
z = X(3);
vx = X(4);
vy = X(5);
vz = X(6);
x_in = (x + mu) * c - y * s;
y_in = (x + mu) * s + y * c;
vx_in = (vx - y) * c - (vy + x + mu) * s;
vy_in = (vx - y) * s + (vy + x + mu) * c;
X_in = [x_in; y_in; z; vx_in; vy_in; vz];
end

function vinf_target = local_active_vinf_target_symbolic(y, t_final, problem)
%LOCAL_ACTIVE_VINF_TARGET_SYMBOLIC Build the symbolic ECI v-infinity target for the active target frame.

% Keep the historical fixed-ECI target in the default mode.
if ~isfield(problem, 'target_frame') || ~strcmp(problem.target_frame, 'earth_sun_rth')
    if strcmp(problem.mode, 'reach')
        vinf_target = problem.vinf_dir(:) * y(problem.vinf_mag_idx);
    else
        vinf_target = problem.vinf_target(:);
    end
    return
end

% Read fixed or free thetaSun0 from the problem payload.
theta_sun0 = problem.theta_sun0_value;
if isfield(problem, 'theta_sun0_idx') && ~isempty(problem.theta_sun0_idx)
    theta_sun0 = y(problem.theta_sun0_idx);
end

% Build the fixed RTH vector or the reachability magnitude-scaled vector.
v_rth = problem.vinf_target_rth(:);
if strcmp(problem.mode, 'reach')
    v_rth_norm = sqrt(sum(v_rth.^2));
    v_rth = v_rth / v_rth_norm * y(problem.vinf_mag_idx);
end

% Convert the active Earth-Sun RTH vector into ECI at Earth perigee.
theta_sun_em = theta_sun0 + problem.ws_sun_em_nd * t_final;
theta_sun_eci = theta_sun_em + t_final;
uR = [-cos(theta_sun_eci); -sin(theta_sun_eci); 0];
uT = [ sin(theta_sun_eci); -cos(theta_sun_eci); 0];
uH = [0; 0; 1];
vinf_target = v_rth(1) * uR + v_rth(2) * uT + v_rth(3) * uH;
end

function vinf_out = local_earth_vinf_out_symbolic(rp_vec, vp_vec, muE)
%LOCAL_EARTH_VINF_OUT_SYMBOLIC Build the symbolic outgoing v-infinity vector at Earth periapsis.

rp = sqrt(sumsqr(rp_vec));
vp = sqrt(sumsqr(vp_vec));
energy = 0.5 * vp^2 - muE / rp;
vinf = sqrt(2 * energy);
delta_half = asin(muE / (muE + rp * vinf^2));
e_r = rp_vec / rp;
e_t = vp_vec / vp;
vinf_out = vinf * (-e_r * sin(delta_half) + e_t * cos(delta_half));
end
