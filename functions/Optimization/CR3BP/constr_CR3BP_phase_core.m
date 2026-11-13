function [c, ceq, Jc, Jceq] = constr_CR3BP_phase_core(yy, problem_in)
%CONSTR_CR3BP_PHASE_CORE Evaluate shared multiple-shooting constraints and Jacobians.
%
% PROTOTYPE:
%   [c, ceq, Jc, Jceq] = constr_CR3BP_phase_core(yy, problem_in)
%
% DESCRIPTION:
%   Enforce segment continuity, departure state/burn, Earth and Moon periapsis,
%   outgoing v-infinity, time limits, and optional maneuver budget/epoch laws.
%   Rebuild absolute node times from the active layout. Segment STMs map
%   initial perturbations to the segment end; returned Jacobians are
%   transposed into MATLAB nonlinear-constraint convention.
%
% INPUT:
%   yy          [nvars x 1] Decision vector from RPO2DS_transfer setup.
%               First 9*n entries are node blocks [X(6); U(3)] in DU/VU.
%               Subsequent time variables use TU; optional RPO phase uses TU,
%               Sun phase rad, and reachable v-infinity magnitude km/s.
%   problem_in  problem_data from the transfer setup, or the supported numeric
%               legacy constant vector. Records n, time/phase indices, active
%               impulse mask, model/scaling constants, targets, periapsis node,
%               altitude/time bounds, and opt/reach mode.
%
% OUTPUT:
%   c           [Nc x 1] Inequalities, feasible when c <= 0.
%   ceq         [Ne x 1] Equalities, feasible when ceq = 0.
%   Jc, Jceq    Optional dense [nvars x Nc/Ne] derivative matrices. Row counts
%               depend on the active timing/phase/reachability policies.
%               Constraint components retain their normalized/physical units.

problem = local_normalize_problem(problem_in, numel(yy));

% Read the scalar constants used across all constraints.
mu = problem.mu;
rminEarth = problem.rminEarth;
rminMoon = problem.rminMoon;
muE = problem.muE;
DU = problem.DU;
VU = problem.VU;
n = problem.n;
maxTransferTOF = problem.maxTransferTOF;
maxTotalTime = problem.maxTotalTime;
idx_peri_moon = problem.idx_peri_moon;
I_2b = problem.I_2b(:);
time_idx = problem.time_idx(:);
dt_idx = problem.dt_idx(:);
t1_idx = problem.t1_idx;
has_final_anchor = strcmp(problem.time_variable_mode, 'propagation_times') && strcmp(problem.time_anchor_mode, 'final');
has_final_time_lock = strcmp(problem.time_variable_mode, 'propagation_times') && problem.lock_t_nodeN && ~has_final_anchor;
has_explicit_phase_variable = problem.use_phase_variable && problem.has_phase_variable && strcmp(problem.phase_source, 'variable');
has_time_derived_phase = problem.use_phase_variable && strcmp(problem.phase_source, 'initial_time');
has_operational_timing = problem.use_phase_variable && ...
    (strcmp(problem.operational_mode, 'initial_state') || strcmp(problem.operational_mode, 'fixed_arrival'));
has_departure_epoch_coupling = has_explicit_phase_variable && ...
    (strcmp(problem.operational_mode, 'initial_state') || strcmp(problem.operational_mode, 'fixed_arrival'));
has_departure_after_initial = has_operational_timing && ~has_final_anchor;
has_operational_final_epoch = problem.use_phase_variable && strcmp(problem.operational_mode, 'fixed_arrival') && ~has_final_anchor;
has_operational_total_time = problem.use_phase_variable && strcmp(problem.operational_mode, 'initial_state');
has_initial_wait_cap = problem.use_phase_variable && strcmp(problem.operational_mode, 'initial_state') && ...
    isfield(problem, 'initial_wait_cap') && ~isempty(problem.initial_wait_cap);

% Allocate the constraint vectors with the active row counts.
ceq = zeros(6 * n + 5 + has_final_time_lock + has_departure_epoch_coupling + has_operational_final_epoch, 1);
if strcmp(problem.mode, 'opt')
    c = zeros(n + 4 + has_operational_total_time + has_departure_after_initial + has_initial_wait_cap, 1);
else
    c = zeros(n + 5 + has_operational_total_time + has_departure_after_initial + has_initial_wait_cap, 1);
end

transfer_cap_row = n + 2;
moon_distance_row = n + 3;
moon_stationarity_row = n + 4;
if has_operational_total_time
    total_time_row = n + 5;
else
    total_time_row = [];
end
if has_departure_after_initial
    departure_after_row = n + 5 + has_operational_total_time;
else
    departure_after_row = [];
end
if has_initial_wait_cap
    initial_wait_cap_row = n + 5 + has_operational_total_time + has_departure_after_initial;
else
    initial_wait_cap_row = [];
end
if strcmp(problem.mode, 'reach')
    dv_budget_row = n + 5 + has_operational_total_time + has_departure_after_initial + has_initial_wait_cap;
else
    dv_budget_row = [];
end

% Allocate dense Jacobians before transposing them at the end.
if nargout > 2
    Jceq = zeros(numel(ceq), numel(yy));
    Jc = zeros(numel(c), numel(yy));
else
    Jceq = [];
    Jc = [];
end

% Rebuild the absolute node times from the active packed layout.
[t_nodes, dt_nodes] = local_decode_time_nodes(problem, yy);
transfer_tof = t_nodes(end) - t_nodes(1);
target_info = local_evaluate_vinf_target(problem, yy, t_nodes(end));
vinf_target = target_info.vinf_target(:);

% Enforce either monotone node times or nonnegative propagation durations.
if strcmp(problem.time_variable_mode, 'propagation_times')
    c(1:n-1) = -dt_nodes;
    c(transfer_cap_row) = sum(dt_nodes) - maxTransferTOF;
    if nargout > 2
        for i = 1:n-1
            Jc(i, dt_idx(i)) = -1;
        end
        Jc(transfer_cap_row, dt_idx(:)) = 1;
    end
else
    c(1:n-1) = t_nodes(1:n-1) - t_nodes(2:n);
    c(transfer_cap_row) = t_nodes(end) - t_nodes(1) - maxTransferTOF;
    if nargout > 2
        for i = 1:n-1
            Jc(i, time_idx(i)) = 1;
            Jc(i, time_idx(i + 1)) = -1;
        end
        Jc(transfer_cap_row, time_idx(1)) = -1;
        Jc(transfer_cap_row, time_idx(end)) = 1;
    end
end

% Enforce the total mission-time cap from the onboard initial condition when requested.
if has_operational_total_time
    c(total_time_row) = t_nodes(end) - problem.initial_epoch_time - maxTotalTime;
    if nargout > 2
        if strcmp(problem.time_variable_mode, 'propagation_times')
            Jc(total_time_row, dt_idx(:)) = 1;
            if ~isempty(t1_idx)
                Jc(total_time_row, t1_idx) = 1;
            end
        else
            Jc(total_time_row, time_idx(end)) = 1;
        end
    end
end

% Keep operational departures after the known onboard initial epoch.
if has_departure_after_initial
    c(departure_after_row) = problem.initial_epoch_time - t_nodes(1);
    if nargout > 2
        if strcmp(problem.time_variable_mode, 'propagation_times')
            if ~isempty(t1_idx)
                Jc(departure_after_row, t1_idx) = -1;
            end
        else
            Jc(departure_after_row, time_idx(1)) = -1;
        end
    end
end

% Limit the RPO waiting time to two RPO periods in initial-state mode.
if has_initial_wait_cap
    c(initial_wait_cap_row) = t_nodes(1) - problem.initial_epoch_time - problem.initial_wait_cap;
    if nargout > 2
        if strcmp(problem.time_variable_mode, 'propagation_times')
            if ~isempty(t1_idx)
                Jc(initial_wait_cap_row, t1_idx) = 1;
            end
        else
            Jc(initial_wait_cap_row, time_idx(1)) = 1;
        end
    end
end

% Enforce the delta-v budget only in reachability mode.
if strcmp(problem.mode, 'reach')
    deltavMag = 0;
    for j = 1:n
        base = 9 * (j - 1);
        u = yy(base + 7);
        w = yy(base + 8);
        s = yy(base + 9);
        deltavMag = deltavMag + (u^2 + w^2 + s^2);
        if nargout > 2
            Jc(dv_budget_row, base + 7) = 2 * u;
            Jc(dv_budget_row, base + 8) = 2 * w;
            Jc(dv_budget_row, base + 9) = 2 * s;
        end
    end
    c(dv_budget_row) = deltavMag - problem.dv_budget;
end

% Build the departure boundary state from a fixed, explicit, or time-derived RPO phase.
[x_ini, dxini_dtheta] = local_get_departure_state(problem, yy, mu, t_nodes(1));

% Read the first node and impose the departure boundary constraints.
x1 = yy(1);
y1 = yy(2);
z1 = yy(3);
vx1 = yy(4);
vy1 = yy(5);
vz1 = yy(6);
u1 = yy(7);
w1 = yy(8);
s1 = yy(9);
ceq(6 * n - 5) = x1 - x_ini(1);
ceq(6 * n - 4) = y1 - x_ini(2);
ceq(6 * n - 3) = z1 - x_ini(3);
check1 = ~isempty(intersect(I_2b, 1));
if ~check1
    ceq(6 * n - 2) = (vx1 - x_ini(4)) - (u1^2 - w1^2 - s1^2);
    ceq(6 * n - 1) = (vy1 - x_ini(5)) - (2 * u1 * w1);
    ceq(6 * n) = (vz1 - x_ini(6)) - (2 * u1 * s1);
else
    ceq(6 * n - 2) = (vx1 - x_ini(4)) - (2 * u1 * w1);
    ceq(6 * n - 1) = (vy1 - x_ini(5)) - (u1^2 - w1^2 - s1^2);
    ceq(6 * n) = (vz1 - x_ini(6)) - (2 * u1 * s1);
end
if nargout > 2
    Jceq(6 * n - 5, 1) = 1;
    Jceq(6 * n - 4, 2) = 1;
    Jceq(6 * n - 3, 3) = 1;
    if has_explicit_phase_variable
        Jceq(6 * n - 5 : 6 * n, problem.phase_theta_idx) = -dxini_dtheta;
    elseif has_time_derived_phase
        if strcmp(problem.time_variable_mode, 'propagation_times')
            if has_final_anchor
                Jceq(6 * n - 5 : 6 * n, dt_idx(:)) = repmat(dxini_dtheta, 1, numel(dt_idx));
            elseif ~isempty(t1_idx)
                Jceq(6 * n - 5 : 6 * n, t1_idx) = -dxini_dtheta;
            end
        else
            Jceq(6 * n - 5 : 6 * n, time_idx(1)) = -dxini_dtheta;
        end
    end
    if ~check1
        Jceq(6 * n - 2, 4) = 1;
        Jceq(6 * n - 2, 7) = -2 * u1;
        Jceq(6 * n - 2, 8) = 2 * w1;
        Jceq(6 * n - 2, 9) = 2 * s1;
        Jceq(6 * n - 1, 5) = 1;
        Jceq(6 * n - 1, 7) = -2 * w1;
        Jceq(6 * n - 1, 8) = -2 * u1;
        Jceq(6 * n, 6) = 1;
        Jceq(6 * n, 7) = -2 * s1;
        Jceq(6 * n, 9) = -2 * u1;
    else
        Jceq(6 * n - 2, 4) = 1;
        Jceq(6 * n - 2, 7) = -2 * w1;
        Jceq(6 * n - 2, 8) = -2 * u1;
        Jceq(6 * n - 1, 5) = 1;
        Jceq(6 * n - 1, 7) = -2 * u1;
        Jceq(6 * n - 1, 8) = 2 * w1;
        Jceq(6 * n - 1, 9) = 2 * s1;
        Jceq(6 * n, 6) = 1;
        Jceq(6 * n, 7) = -2 * s1;
        Jceq(6 * n, 9) = -2 * u1;
    end
end

% Enforce the multiple-shooting continuity defects across every segment.
for i = 1:n-1
    check = ~isempty(intersect(I_2b, i + 1));
    xi = yy(9 * (i - 1) + 1);
    yi = yy(9 * (i - 1) + 2);
    zi = yy(9 * (i - 1) + 3);
    vxi = yy(9 * (i - 1) + 4);
    vyi = yy(9 * (i - 1) + 5);
    vzi = yy(9 * (i - 1) + 6);
    Xi = [xi, yi, zi, vxi, vyi, vzi]';
    xip1 = yy(9 * i + 1);
    yip1 = yy(9 * i + 2);
    zip1 = yy(9 * i + 3);
    vxip1 = yy(9 * i + 4);
    vyip1 = yy(9 * i + 5);
    vzip1 = yy(9 * i + 6);
    uip1 = yy(9 * i + 7);
    wip1 = yy(9 * i + 8);
    sip1 = yy(9 * i + 9);
    ti = t_nodes(i);
    tip1 = t_nodes(i + 1);

    % Propagate one segment with or without the STM as needed.
    if nargout < 3
        x_temp = propagate(Xi, ti, tip1, mu, '3D', 'CR3BP', 'false');
    else
        [x_temp, PHI] = propagate(Xi, ti, tip1, mu, '3D', 'CR3BP', 'true');
    end

    % Enforce the position and velocity continuity defects.
    ceq(3 * (i - 1) + 1) = x_temp(1) - xip1;
    ceq(3 * (i - 1) + 2) = x_temp(2) - yip1;
    ceq(3 * (i - 1) + 3) = x_temp(3) - zip1;
    if ~check
        ceq(3 * n - 5 + 3 * i) = (vxip1 - x_temp(4)) - (uip1^2 - wip1^2 - sip1^2);
        ceq(3 * n - 4 + 3 * i) = (vyip1 - x_temp(5)) - (2 * uip1 * wip1);
        ceq(3 * n - 3 + 3 * i) = (vzip1 - x_temp(6)) - (2 * uip1 * sip1);
    else
        ceq(3 * n - 5 + 3 * i) = (vxip1 - x_temp(4)) - (2 * uip1 * wip1);
        ceq(3 * n - 4 + 3 * i) = (vyip1 - x_temp(5)) - (uip1^2 - wip1^2 - sip1^2);
        ceq(3 * n - 3 + 3 * i) = (vzip1 - x_temp(6)) - (2 * uip1 * sip1);
    end

    % Build the analytic Jacobian for the local segment defect rows.
    if nargout > 2
        fi = rhs_CR3BP_3D(ti, Xi, mu);
        f_temp = rhs_CR3BP_3D(tip1, x_temp, mu);
        pos_rows = 3 * (i - 1) + 1 : 3 * i;
        vel_rows = 3 * n - 5 + 3 * i : 3 * n - 3 + 3 * i;
        Jceq(pos_rows, 9 * (i - 1) + 1 : 9 * (i - 1) + 6) = PHI(1:3, 1:6);
        Jceq(pos_rows, 9 * i + 1 : 9 * i + 3) = -eye(3);
        Jceq(vel_rows, 9 * i + 4 : 9 * i + 6) = eye(3);
        Jceq(vel_rows, 9 * (i - 1) + 1 : 9 * (i - 1) + 6) = -PHI(4:6, 1:6);
        if strcmp(problem.time_variable_mode, 'propagation_times')
            Jceq(pos_rows, dt_idx(i)) = f_temp(1:3, 1);
            Jceq(vel_rows, dt_idx(i)) = -f_temp(4:6, 1);
        else
            Jceq(pos_rows, time_idx(i)) = -PHI(1:3, 1:6) * fi(1:6, 1);
            Jceq(pos_rows, time_idx(i + 1)) = f_temp(1:3, 1);
            Jceq(vel_rows, time_idx(i)) = PHI(4:6, 1:6) * fi(1:6, 1);
            Jceq(vel_rows, time_idx(i + 1)) = -f_temp(4:6, 1);
        end
        if ~check
            Jceq(vel_rows, 9 * i + 7 : 9 * i + 9) = 2 * [-uip1, wip1, sip1; -wip1, -uip1, 0; -sip1, 0, -uip1];
        else
            Jceq(vel_rows, 9 * i + 7 : 9 * i + 9) = 2 * [-wip1, -uip1, 0; -uip1, wip1, sip1; -sip1, 0, -uip1];
        end
    end
end

% Enforce the lunar-periapsis minimum-distance and stationarity conditions.
X_fb = yy(9 * idx_peri_moon - 8 : 9 * idx_peri_moon - 3)';
[r2, r2dot, r2ddot] = get_r2ddot(X_fb, mu);
c(moon_distance_row) = rminMoon^2 - r2^2;
ceq(6 * n + 5) = r2dot;
c(moon_stationarity_row) = -r2ddot;
if nargout > 2
    event_constants = local_build_event_constants(problem, vinf_target);
    [dr2dot_dt, dr2dot_dx, dr2dot_dy, dr2dot_dz, ...
        dr2dot_dvx, dr2dot_dvy, dr2dot_dvz, ...
        dr2ddot_dt, dr2ddot_dx, dr2ddot_dy, dr2ddot_dz, ...
        dr2ddot_dvx, dr2ddot_dvy, dr2ddot_dvz] = get_dr2ddot_dtX(X_fb, event_constants);
    Jc(moon_distance_row, 9 * idx_peri_moon - 8) = -2 * (X_fb(1) - 1 + mu);
    Jc(moon_distance_row, 9 * idx_peri_moon - 7) = -2 * X_fb(2);
    Jc(moon_distance_row, 9 * idx_peri_moon - 6) = -2 * X_fb(3);
    Jceq(6 * n + 5, 9 * idx_peri_moon - 8) = dr2dot_dx;
    Jceq(6 * n + 5, 9 * idx_peri_moon - 7) = dr2dot_dy;
    Jceq(6 * n + 5, 9 * idx_peri_moon - 6) = dr2dot_dz;
    Jceq(6 * n + 5, 9 * idx_peri_moon - 5) = dr2dot_dvx;
    Jceq(6 * n + 5, 9 * idx_peri_moon - 4) = dr2dot_dvy;
    Jceq(6 * n + 5, 9 * idx_peri_moon - 3) = dr2dot_dvz;
    Jc(moon_stationarity_row, 9 * idx_peri_moon - 8) = -dr2ddot_dx;
    Jc(moon_stationarity_row, 9 * idx_peri_moon - 7) = -dr2ddot_dy;
    Jc(moon_stationarity_row, 9 * idx_peri_moon - 6) = -dr2ddot_dz;
    Jc(moon_stationarity_row, 9 * idx_peri_moon - 5) = -dr2ddot_dvx;
    Jc(moon_stationarity_row, 9 * idx_peri_moon - 4) = -dr2ddot_dvy;
    Jc(moon_stationarity_row, 9 * idx_peri_moon - 3) = -dr2ddot_dvz;
    if strcmp(problem.time_variable_mode, 'propagation_times')
        if has_final_anchor
            if idx_peri_moon <= numel(dt_idx)
                Jceq(6 * n + 5, dt_idx(idx_peri_moon:end)) = -dr2dot_dt;
                Jc(moon_stationarity_row, dt_idx(idx_peri_moon:end)) = dr2ddot_dt;
            end
        elseif idx_peri_moon > 1
            Jceq(6 * n + 5, dt_idx(1:idx_peri_moon-1)) = dr2dot_dt;
            Jc(moon_stationarity_row, dt_idx(1:idx_peri_moon-1)) = -dr2ddot_dt;
        end
        if ~has_final_anchor && ~isempty(t1_idx)
            Jceq(6 * n + 5, t1_idx) = dr2dot_dt;
            Jc(moon_stationarity_row, t1_idx) = -dr2ddot_dt;
        end
    else
        Jceq(6 * n + 5, time_idx(idx_peri_moon)) = dr2dot_dt;
        Jc(moon_stationarity_row, time_idx(idx_peri_moon)) = -dr2ddot_dt;
    end
end

% Enforce the Earth-periapsis boundary and final vInf matching.
Xn = yy(9 * n - 8 : 9 * n - 3)';
tn = t_nodes(end);
[r1, r1dot, r1ddot] = get_r1ddot(Xn, mu);
c(n + 1) = rminEarth^2 - r1^2;
ceq(6 * n + 1) = r1dot;
c(n) = -r1ddot;
Xn_in = rot2in(Xn, tn, mu);
Xn_in(1:3) = Xn_in(1:3) * DU;
Xn_in(4:6) = Xn_in(4:6) * VU;
[~, ~, vinf_out_DVaft] = earth_flyby_for_optimization(Xn_in(1:3), Xn_in(4:6), muE);
ceq(6 * n + 2) = vinf_out_DVaft(1) - vinf_target(1);
ceq(6 * n + 3) = vinf_out_DVaft(2) - vinf_target(2);
ceq(6 * n + 4) = vinf_out_DVaft(3) - vinf_target(3);
if nargout > 2
    event_constants = local_build_event_constants(problem, vinf_target);
    [dr1dot_dt, dr1dot_dx, dr1dot_dy, dr1dot_dz, dr1dot_dvx, dr1dot_dvy, dr1dot_dvz, ...
        dr1ddot_dt, dr1ddot_dx, dr1ddot_dy, dr1ddot_dz, dr1ddot_dvx, dr1ddot_dvy, dr1ddot_dvz] = get_dr1ddot_dtX(Xn, event_constants);
    Jc(n + 1, 9 * n - 8) = -2 * (Xn(1) + mu);
    Jc(n + 1, 9 * n - 7) = -2 * Xn(2);
    Jc(n + 1, 9 * n - 6) = -2 * Xn(3);
    Jceq(6 * n + 1, 9 * n - 8) = dr1dot_dx;
    Jceq(6 * n + 1, 9 * n - 7) = dr1dot_dy;
    Jceq(6 * n + 1, 9 * n - 6) = dr1dot_dz;
    Jceq(6 * n + 1, 9 * n - 5) = dr1dot_dvx;
    Jceq(6 * n + 1, 9 * n - 4) = dr1dot_dvy;
    Jceq(6 * n + 1, 9 * n - 3) = dr1dot_dvz;
    Jc(n, 9 * n - 8) = -dr1ddot_dx;
    Jc(n, 9 * n - 7) = -dr1ddot_dy;
    Jc(n, 9 * n - 6) = -dr1ddot_dz;
    Jc(n, 9 * n - 5) = -dr1ddot_dvx;
    Jc(n, 9 * n - 4) = -dr1ddot_dvy;
    Jc(n, 9 * n - 3) = -dr1ddot_dvz;
    dvinfdXn_block = local_dvinfdXn(Xn, tn, mu, muE, DU, VU);
    Jceq(6 * n + 2 : 6 * n + 4, 9 * n - 8 : 9 * n - 3) = dvinfdXn_block(:, 1:6);
    if strcmp(problem.time_variable_mode, 'propagation_times')
        if ~has_final_anchor
            Jceq(6 * n + 1, dt_idx(:)) = dr1dot_dt;
            Jc(n, dt_idx(:)) = -dr1ddot_dt;
            Jceq(6 * n + 2 : 6 * n + 4, dt_idx(:)) = repmat(dvinfdXn_block(:, 7) - target_info.dtarget_dt, 1, numel(dt_idx));
        end
        if ~has_final_anchor && ~isempty(t1_idx)
            Jceq(6 * n + 1, t1_idx) = dr1dot_dt;
            Jc(n, t1_idx) = -dr1ddot_dt;
            Jceq(6 * n + 2 : 6 * n + 4, t1_idx) = dvinfdXn_block(:, 7) - target_info.dtarget_dt;
        end
    else
        Jceq(6 * n + 1, time_idx(end)) = dr1dot_dt;
        Jc(n, time_idx(end)) = -dr1ddot_dt;
        Jceq(6 * n + 2 : 6 * n + 4, time_idx(end)) = dvinfdXn_block(:, 7) - target_info.dtarget_dt;
    end
    if isfield(problem, 'theta_sun0_idx') && ~isempty(problem.theta_sun0_idx)
        Jceq(6 * n + 2 : 6 * n + 4, problem.theta_sun0_idx) = -target_info.dtarget_dtheta0;
    end
    if strcmp(problem.mode, 'reach')
        Jceq(6 * n + 2 : 6 * n + 4, problem.vinf_mag_idx) = -target_info.dtarget_dmag;
    end
end

% Enforce the prescribed final epoch when propagation_times does not store t(N) directly.
if has_final_time_lock
    final_time_row = 6 * n + 6;
    ceq(final_time_row) = t_nodes(end) - problem.t_nodeN_fixed;
    if nargout > 2
        Jceq(final_time_row, dt_idx(:)) = 1;
        if ~isempty(t1_idx)
            Jceq(final_time_row, t1_idx) = 1;
        end
    end
end

% Enforce the absolute departure epoch implied by the onboard initial state.
if has_departure_epoch_coupling
    coupling_row = 6 * n + 5 + has_final_time_lock + 1;
    ceq(coupling_row) = t_nodes(1) - problem.initial_epoch_time - ...
        yy(problem.phase_theta_idx) + problem.initial_rpo_phase;
    if nargout > 2
        Jceq(coupling_row, problem.phase_theta_idx) = -1;
        if strcmp(problem.time_variable_mode, 'propagation_times')
            if ~isempty(t1_idx)
                Jceq(coupling_row, t1_idx) = 1;
            end
        else
            Jceq(coupling_row, time_idx(1)) = 1;
        end
    end
end

% Enforce the fixed absolute final epoch in the fixed-arrival operational mode.
if has_operational_final_epoch
    final_epoch_row = 6 * n + 5 + has_final_time_lock + has_departure_epoch_coupling + 1;
    ceq(final_epoch_row) = t_nodes(end) - problem.final_epoch_time;
    if nargout > 2
        if strcmp(problem.time_variable_mode, 'propagation_times')
            Jceq(final_epoch_row, dt_idx(:)) = 1;
            if ~isempty(t1_idx)
                Jceq(final_epoch_row, t1_idx) = 1;
            end
        else
            Jceq(final_epoch_row, time_idx(end)) = 1;
        end
    end
end

% Return the Jacobians in the orientation expected by fmincon.
if nargout > 2
    Jceq = Jceq';
    Jc = Jc';
end
end

function problem = local_normalize_problem(problem_in, nvars)
%LOCAL_NORMALIZE_PROBLEM Normalize numeric legacy inputs and new structs into one shared layout.

if isstruct(problem_in)
    problem = problem_in;
    problem.time_variable_mode = char(lower(string(problem.time_variable_mode)));
    if ~isfield(problem, 'time_anchor_mode') || isempty(problem.time_anchor_mode)
        problem.time_anchor_mode = 'initial';
    else
        problem.time_anchor_mode = char(lower(string(problem.time_anchor_mode)));
    end
    if ~isfield(problem, 'operational_mode') || isempty(problem.operational_mode)
        problem.operational_mode = 'standard';
    else
        problem.operational_mode = char(lower(string(problem.operational_mode)));
    end
    if ~isfield(problem, 'geometry_mode') || isempty(problem.geometry_mode)
        problem.geometry_mode = 'free_geometry';
    end
    if ~isfield(problem, 'time_idx') || isempty(problem.time_idx)
        if strcmp(problem.time_variable_mode, 'propagation_times')
            problem.time_idx = (9 * problem.n + 1 : 10 * problem.n - 1).';
        else
            problem.time_idx = (9 * problem.n + 1 : 10 * problem.n).';
        end
    end
    if ~isfield(problem, 'dt_idx') || isempty(problem.dt_idx)
        if strcmp(problem.time_variable_mode, 'propagation_times')
            problem.dt_idx = problem.time_idx(:);
        else
            problem.dt_idx = [];
        end
    end
    if ~isfield(problem, 't1_idx') || isempty(problem.t1_idx)
        problem.t1_idx = [];
    end
    if ~isfield(problem, 'final_time_idx') || isempty(problem.final_time_idx)
        problem.final_time_idx = problem.time_idx(end);
    end
    if ~isfield(problem, 't0_anchor') || isempty(problem.t0_anchor)
        problem.t0_anchor = [];
    end
    if ~isfield(problem, 'lock_t_node1') || isempty(problem.lock_t_node1)
        problem.lock_t_node1 = false;
    end
    if ~isfield(problem, 'lock_t_nodeN') || isempty(problem.lock_t_nodeN)
        problem.lock_t_nodeN = false;
    end
    if ~isfield(problem, 't_node1_fixed')
        problem.t_node1_fixed = [];
    end
    if ~isfield(problem, 't_nodeN_fixed')
        problem.t_nodeN_fixed = [];
    end
    if ~isfield(problem, 'has_t1_variable') || isempty(problem.has_t1_variable)
        problem.has_t1_variable = ~isempty(problem.t1_idx);
    end
    if ~isfield(problem, 'use_phase_variable') || isempty(problem.use_phase_variable)
        problem.use_phase_variable = isfield(problem, 'phase_theta_idx') && ~isempty(problem.phase_theta_idx);
    end
    if ~isfield(problem, 'phase_theta_idx')
        problem.phase_theta_idx = [];
    end
    if ~isfield(problem, 'has_phase_variable') || isempty(problem.has_phase_variable)
        problem.has_phase_variable = ~isempty(problem.phase_theta_idx);
    end
    if ~isfield(problem, 'phase_source') || isempty(problem.phase_source)
        if problem.has_phase_variable
            problem.phase_source = 'variable';
        else
            problem.phase_source = 'fixed';
        end
    else
        problem.phase_source = char(lower(string(problem.phase_source)));
    end
    if ~isfield(problem, 'nvars') || isempty(problem.nvars)
        problem.nvars = nvars;
    end
    if ~isfield(problem, 'phase_origin') || isempty(problem.phase_origin)
        problem.phase_origin = 0;
    end
    if ~isfield(problem, 'phase_period') || isempty(problem.phase_period)
        problem.phase_period = [];
    end
    if ~isfield(problem, 'initial_rpo_phase') || isempty(problem.initial_rpo_phase)
        problem.initial_rpo_phase = 0;
    end
    if ~isfield(problem, 'initial_epoch_time') || isempty(problem.initial_epoch_time)
        problem.initial_epoch_time = 0;
    end
    if ~isfield(problem, 'final_epoch_time')
        problem.final_epoch_time = [];
    end
    if ~isfield(problem, 'available_time_nd')
        problem.available_time_nd = [];
    end
    if ~isfield(problem, 'maxTransferTOF') || isempty(problem.maxTransferTOF)
        problem.maxTransferTOF = problem.maxTOF;
    end
    if ~isfield(problem, 'maxTotalTime')
        problem.maxTotalTime = [];
    end
    if ~isfield(problem, 'target_frame') || isempty(problem.target_frame)
        problem.target_frame = 'eci';
    else
        problem.target_frame = char(lower(string(problem.target_frame)));
    end
    if ~isfield(problem, 'theta_sun0_mode') || isempty(problem.theta_sun0_mode)
        problem.theta_sun0_mode = 'specified';
    end
    if ~isfield(problem, 'theta_sun0_value') || isempty(problem.theta_sun0_value)
        problem.theta_sun0_value = 0;
    end
    if ~isfield(problem, 'theta_sun0_idx')
        problem.theta_sun0_idx = [];
    end
    if ~isfield(problem, 'ws_sun_em_nd') || isempty(problem.ws_sun_em_nd)
        problem.ws_sun_em_nd = -9.25195985e-1;
    end
    if ~isfield(problem, 'vinf_target_rth') || isempty(problem.vinf_target_rth)
        problem.vinf_target_rth = [0; 3; 0];
    end
    if ~isfield(problem, 'vinf_mag_idx') || isempty(problem.vinf_mag_idx)
        problem.vinf_mag_idx = nvars;
    end
    return
end

% Rebuild the old numeric payload as the legacy node-time formulation.
constants = problem_in(:);
n = constants(16);
problem = struct();
problem.mu = constants(1);
problem.rminEarth = constants(2);
problem.rminMoon = constants(3);
problem.muE = constants(4);
problem.DU = constants(14);
problem.VU = constants(15);
problem.TU = [];
problem.n = n;
problem.maxTOF = constants(17);
problem.idx_peri_moon = constants(18);
problem.x_ini_fixed = constants(8:13);
problem.x0_rpo_base = constants(8:13);
problem.use_phase_variable = false;
problem.has_phase_variable = false;
problem.phase_source = 'fixed';
problem.phase_theta_idx = [];
problem.phase_origin = 0;
problem.phase_period = [];
problem.time_variable_mode = 'node_times';
problem.time_anchor_mode = 'initial';
problem.geometry_mode = 'free_geometry';
problem.t0_anchor = [];
problem.time_idx = (9 * n + 1 : 10 * n).';
problem.dt_idx = [];
problem.t1_idx = [];
problem.time_count = n;
problem.final_time_idx = 10 * n;
problem.lock_t_node1 = false;
problem.lock_t_nodeN = false;
problem.t_node1_fixed = [];
problem.t_nodeN_fixed = [];
problem.has_t1_variable = false;
problem.operational_mode = 'standard';
problem.initial_rpo_phase = 0;
problem.initial_epoch_time = 0;
problem.final_epoch_time = [];
problem.available_time_nd = [];
problem.maxTransferTOF = problem.maxTOF;
problem.maxTotalTime = [];
problem.nvars = nvars;
if nvars == 10 * n + 1
    problem.mode = 'reach';
    problem.vinf_dir = constants(5:7);
    problem.dv_budget = constants(19);
    if length(constants) > 19
        problem.I_2b = constants(20:end);
    else
        problem.I_2b = [];
    end
else
    problem.mode = 'opt';
    problem.vinf_target = constants(5:7);
    if length(constants) > 18
        problem.I_2b = constants(19:end);
    else
        problem.I_2b = [];
    end
end
end

function [t_nodes, dt_nodes] = local_decode_time_nodes(problem, yy)
%LOCAL_DECODE_TIME_NODES Rebuild the absolute times and segment durations from the active layout.

if strcmp(problem.time_variable_mode, 'propagation_times')
    dt_nodes = yy(problem.dt_idx(:));
    if strcmp(problem.time_anchor_mode, 'final')
        t_nodes = problem.t0_anchor - [flipud(cumsum(flipud(dt_nodes))); 0];
    elseif ~isempty(problem.t1_idx)
        t1 = yy(problem.t1_idx);
        t_nodes = t1 + [0; cumsum(dt_nodes)];
    else
        t1 = problem.t0_anchor;
        t_nodes = t1 + [0; cumsum(dt_nodes)];
    end
else
    t_nodes = yy(problem.time_idx(:));
    dt_nodes = diff(t_nodes);
end
end

function [x_ini, dxini_dtheta] = local_get_departure_state(problem, yy, mu, t1)
%LOCAL_GET_DEPARTURE_STATE Build the departure boundary state for the fixed-phase or phase-aware case.

if problem.use_phase_variable && strcmp(problem.phase_source, 'initial_time')
    theta_raw = problem.initial_rpo_phase + t1 - problem.initial_epoch_time;
    theta_rpo = local_wrap_phase_value(theta_raw, problem.phase_origin, problem.phase_period);
    x_ini = propagate(problem.x0_rpo_base, 0, theta_rpo, mu, '3D', 'CR3BP', 'false');
    % Phase is elapsed time: use the autonomous flow, not differences across
    % the period wrap where the numerical orbit has a nonzero closure error.
    dxini_dtheta = rhs_CR3BP_3D(theta_rpo, x_ini, mu);
elseif problem.use_phase_variable && problem.has_phase_variable
    theta_raw = yy(problem.phase_theta_idx);
    theta_rpo = local_wrap_phase_value(theta_raw, problem.phase_origin, problem.phase_period);
    x_ini = propagate(problem.x0_rpo_base, 0, theta_rpo, mu, '3D', 'CR3BP', 'false');
    dxini_dtheta = rhs_CR3BP_3D(theta_rpo, x_ini, mu);
else
    x_ini = problem.x_ini_fixed(:);
    dxini_dtheta = zeros(6, 1);
end
end

function phase_wrapped = local_wrap_phase_value(phase_raw, phase_origin, phase_period)
%LOCAL_WRAP_PHASE_VALUE Wrap one RPO phase onto the database geometry interval.

% Keep fixed-phase cases unchanged when no period metadata is available.
phase_wrapped = phase_raw;
if isempty(phase_period) || ~isfinite(phase_period) || phase_period <= 0
    return
end

% Remove whole-period waits while preserving the phase origin.
phase_wrapped = phase_origin + mod(phase_raw - phase_origin, phase_period);
end

function event_constants = local_build_event_constants(problem, vinf_target)
%LOCAL_BUILD_EVENT_CONSTANTS Build the legacy event-derivative constant vector used by helper functions.

if problem.use_phase_variable
    x_seed = problem.x0_rpo_base(:);
else
    x_seed = problem.x_ini_fixed(:);
end
event_constants = [problem.mu; problem.rminEarth; problem.rminMoon; problem.muE; ...
    vinf_target(:); x_seed(:); problem.DU; problem.VU; problem.n; ...
    problem.maxTOF; problem.idx_peri_moon; problem.I_2b(:)];
end

function target_info = local_evaluate_vinf_target(problem, yy, t_final)
%LOCAL_EVALUATE_VINF_TARGET Rebuild the active ECI v-infinity target and its key derivatives.

% Use the historical fixed-ECI target unless the user requested RTH.
if ~isfield(problem, 'target_frame') || ~strcmp(problem.target_frame, 'earth_sun_rth')
    if strcmp(problem.mode, 'opt')
        vinf_target = problem.vinf_target(:);
        dtarget_dmag = zeros(3, 1);
    else
        vinf_target = problem.vinf_dir(:) * yy(problem.vinf_mag_idx);
        dtarget_dmag = problem.vinf_dir(:);
    end
    target_info = struct('vinf_target', vinf_target, ...
        'dtarget_dt', zeros(3, 1), ...
        'dtarget_dtheta0', zeros(3, 1), ...
        'dtarget_dmag', dtarget_dmag);
    return
end

% Read fixed or free thetaSun0 from the active problem definition.
theta_sun0 = problem.theta_sun0_value;
if isfield(problem, 'theta_sun0_idx') && ~isempty(problem.theta_sun0_idx)
    theta_sun0 = yy(problem.theta_sun0_idx);
end

% Build the target RTH vector for fixed-target or reachability mode.
v_rth = problem.vinf_target_rth(:);
dtarget_dmag = zeros(3, 1);
if strcmp(problem.mode, 'reach')
    rth_norm = max(norm(v_rth), eps);
    v_rth_dir = v_rth / rth_norm;
    v_rth = v_rth_dir * yy(problem.vinf_mag_idx);
else
    v_rth_dir = zeros(3, 1);
end

% Convert the active Sun geometry into the ECI target vector.
theta_sun_em = theta_sun0 + problem.ws_sun_em_nd * t_final;
[vinf_target, dv_dtheta_eci] = earth_sun_rth_to_eci(v_rth, theta_sun_em, t_final);
if strcmp(problem.mode, 'reach')
    [dtarget_dmag, ~] = earth_sun_rth_to_eci(v_rth_dir, theta_sun_em, t_final);
end

% Pack derivatives with respect to final epoch and thetaSun0.
target_info = struct('vinf_target', vinf_target, ...
    'dtarget_dt', dv_dtheta_eci * (1 + problem.ws_sun_em_nd), ...
    'dtarget_dtheta0', dv_dtheta_eci, ...
    'dtarget_dmag', dtarget_dmag);
end

function dvinfdXn_block = local_dvinfdXn(Xn, tn, mu, muE, DU, VU)
%LOCAL_DVINFDXN Differentiate the outgoing Earth-flyby vInf by finite differences.

% Allocate the state/time Jacobian and use one shared finite-difference step.
dvinfdXn_block = zeros(3, 7);
delta = sqrt(eps);

% Differentiate with respect to the six rotating-frame state components.
for i = 1:6
    Xn_plus = Xn;
    Xn_minus = Xn;
    Xn_plus(i) = Xn_plus(i) + delta;
    Xn_minus(i) = Xn_minus(i) - delta;
    Xin_plus = rot2in(Xn_plus, tn, mu);
    Xin_minus = rot2in(Xn_minus, tn, mu);
    Xin_plus(1:3) = Xin_plus(1:3) * DU;
    Xin_plus(4:6) = Xin_plus(4:6) * VU;
    Xin_minus(1:3) = Xin_minus(1:3) * DU;
    Xin_minus(4:6) = Xin_minus(4:6) * VU;
    [~, ~, vinf_plus] = earth_flyby_for_optimization(Xin_plus(1:3), Xin_plus(4:6), muE);
    [~, ~, vinf_minus] = earth_flyby_for_optimization(Xin_minus(1:3), Xin_minus(4:6), muE);
    dvinfdXn_block(:, i) = (vinf_plus(:) - vinf_minus(:)) / (2 * delta);
end

% Differentiate with respect to the final time separately.
Xin_plus = rot2in(Xn, tn + delta, mu);
Xin_minus = rot2in(Xn, tn - delta, mu);
Xin_plus(1:3) = Xin_plus(1:3) * DU;
Xin_plus(4:6) = Xin_plus(4:6) * VU;
Xin_minus(1:3) = Xin_minus(1:3) * DU;
Xin_minus(4:6) = Xin_minus(4:6) * VU;
[~, ~, vinf_plus] = earth_flyby_for_optimization(Xin_plus(1:3), Xin_plus(4:6), muE);
[~, ~, vinf_minus] = earth_flyby_for_optimization(Xin_minus(1:3), Xin_minus(4:6), muE);
dvinfdXn_block(:, 7) = (vinf_plus(:) - vinf_minus(:)) / (2 * delta);
end
