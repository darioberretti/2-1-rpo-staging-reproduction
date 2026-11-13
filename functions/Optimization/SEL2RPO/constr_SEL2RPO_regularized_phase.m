function [c, ceq, Jc, Jceq] = constr_SEL2RPO_regularized_phase(y, nlp)
%CONSTR_SEL2RPO_REGULARIZED_PHASE Nonlinear constraints for SEL2-to-RPO transfers.
%
% PROTOTYPE:
%   [c, ceq, Jc, Jceq] = constr_SEL2RPO_regularized_phase(y, nlp)
%
% DESCRIPTION:
%   Evaluate BCR4BP shooting defects and periapsis/time constraints using
%   regularized impulses. Boundary-phase sensitivities use the exact CR3BP
%   flow and the SE-to-EM rotating-frame chain rule.
%
% INPUT:
%   y    [10*N-11 x 1] Decision vector: U1, interior (Xj,Uj), UN,
%        segment durations, SEL2 phase and RPO phase. Uj = [u,w,s] defines
%        an Oshima regularized impulse. States are [x,y,z,vx,vy,vz] in
%        EM normalized rotating coordinates; durations and RPO phase use
%        EM time units, and SEL2 phase uses SE time units. [-]
%   nlp  Problem structure from data/inputs/run_04_reference_inputs.mat:
%        node/variable maps, halo and RPO seeds, constants, epoch, integration
%        tolerances and constraint limits. N is the number of shooting nodes.
%
% OUTPUT:
%   c, ceq      Column vectors with c <= 0 and ceq = 0. [-]
%   Jc, Jceq    Sparse [nvars x nc] and [nvars x neq] Jacobians,
%                 variables in rows (fmincon convention). If E is 1 for
%                 an enabled Earth-periapsis node and 0 otherwise,
%                 nc=N+2+2*E and neq=6*(N-1)+1+E.
%
%   Inequalities enforce minimum segment duration, maximum flight time,
%   Moon radius/positive radial acceleration, and optionally Earth
%   radius/positive radial acceleration. Equalities enforce six shooting
%   defects per arc and zero radial speed at the selected periapsides.
%   Boundary sensitivities use exact flow/frame derivatives; radial
%   speed/acceleration gradients retain finite differences. The wrapped
%   halo waiting time has only a local one-sided derivative at a wrap.
%   nlp.I_case2b selects impulses whose squared component lies along y
%   instead of x. nlp.se_secondary_center defaults to earth; Run 04
%   explicitly selects emb to retain its historical frame convention.

% Core dimensions and strict options
N = round(nlp.N);
nvars = numel(y);

% Keep Earth-peri feature as the only on/off branch (valid node index only).
use_earth_peri = ~isempty(nlp.idx_earth_peri_node) && isfinite(nlp.idx_earth_peri_node) && ...
    (nlp.idx_earth_peri_node >= 1) && (nlp.idx_earth_peri_node <= N);

% Decode decision vector
U_nodes = zeros(N, 3);
X_nodes_interior = NaN(N, 6);

% Extract U on all nodes and X on interior nodes.
for j = 1:N
    idxU = nlp.map.node(j).U;
    U_nodes(j, :) = y(idxU(:)).';

    idxX = nlp.map.node(j).X;
    if ~isempty(idxX)
        X_nodes_interior(j, :) = y(idxX(:)).';
    end
end

% Extract segment durations and phase variables.
dt_nodes = y(nlp.map.dt_nodes(:));
theta_sel2 = y(nlp.map.theta_sel2);
theta_rpo = y(nlp.map.theta_rpo);

% Reconstruct timeline from phase and durations
% Convert selected SEL2 phase into wrapped waiting time.
t_start_phase_se_nd = mod(theta_sel2 - nlp.theta_halo0_nd, nlp.sel2_period_nd);
k_se2em = nlp.omega_em_rad_s / nlp.omega_se_rad_s;
t_start_phase_em_nd = k_se2em * t_start_phase_se_nd;

% Differentiate the local waiting-time branch. At a modulo discontinuity
% this is the right-sided derivative; a two-sided derivative does not exist.
dt_start_em_dtheta = k_se2em;

% Build transfer-time node grid from configured origin and dt sequence.
t_nodes = t_start_phase_em_nd + [0; cumsum(dt_nodes)];
transfer_time_nd = sum(dt_nodes);
total_time_nd = t_start_phase_em_nd + transfer_time_nd; %#ok<NASGU>

% Evaluate SEL2 boundary state and sensitivities
% Propagate the halo once to the selected SE phase.
sel2_def = nlp.sel2_traj_phase;
ode_opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
if theta_sel2 == nlp.theta_halo0_nd
    x_se = sel2_def.x0_se_nd(:);
else
    ode_sol_struct = ode113(@(t, xx) rhs_SE3BP_3D(t, xx, nlp.mu_se_nd), ...
        [nlp.theta_halo0_nd, theta_sel2], sel2_def.x0_se_nd(:), ode_opts);
    x_se = ode_sol_struct.y(:, end);
end

% Run 04 explicitly reproduces the March barycenter convention. Other callers
% retain the corrected physical-Earth default.
secondary_center = 'earth';
if isfield(nlp, 'se_secondary_center')
    secondary_center = nlp.se_secondary_center;
end

% Convert the nominal state and differentiate its linear state mapping.
theta_sun_abs = nlp.theta_sun0_rad + nlp.ws_sun_em_nd * t_start_phase_em_nd;
[X_sel2_pre, frame] = se_cr3bp_to_emb_bcr4bp_state(x_se.', theta_sun_abs, ...
    'alpha_plus_theta_sun_const_rad', nlp.alpha_plus_theta_const, ...
    'l_se_km', nlp.l_se_km, 'l_em_km', nlp.l_em_km, ...
    'omega_se_rad_s', nlp.omega_se_rad_s, ...
    'omega_em_rad_s', nlp.omega_em_rad_s, ...
    'mu_se_nd', nlp.mu_se_nd, 'mu_em_nd', nlp.mu_em, ...
    'se_secondary_center', secondary_center);

length_ratio = frame.constants.l_ratio_se_to_em;
rate_ratio = frame.constants.omega_ratio_se_to_em;
alpha = frame.alpha_rad;
rotation = length_ratio * [cos(alpha), sin(alpha); -sin(alpha), cos(alpha)];
quarter_turn = [0, -1; 1, 0];
flow_se = rhs_SE3BP_3D(theta_sel2, x_se, nlp.mu_se_nd);
dposition_dtheta = rotation * flow_se(1:2);
dvelocity_dtheta = rate_ratio * rotation * flow_se(4:5) ...
    - (1 - rate_ratio) * quarter_turn * dposition_dtheta;
dXsel2_dtheta = [dposition_dtheta; length_ratio * flow_se(3); ...
    dvelocity_dtheta; length_ratio * rate_ratio * flow_se(6)];

% At fixed halo phase, only planar position/velocity rotate with EM time.
secondary_x = 0;
if strcmpi(secondary_center, 'earth')
    secondary_x = -nlp.mu_em;
end
dXsel2_dt = nlp.ws_sun_em_nd * [-X_sel2_pre(2); X_sel2_pre(1) - secondary_x; ...
    0; -X_sel2_pre(5); X_sel2_pre(4); 0];

% Evaluate RPO boundary state and sensitivities
% Build one 6D RPO seed state from planar definition.
x0p_rpo = nlp.rpo_traj_phase.x0_rpo_planar(:);
x0_rpo = [x0p_rpo(1); x0p_rpo(2); 0; x0p_rpo(3); x0p_rpo(4); 0];

% The derivative with respect to autonomous orbit phase is its dynamics.
X_rpo_post = propagate(x0_rpo, 0, theta_rpo, nlp.mu_em, '3D', 'CR3BP', 'false');
dXrpo_dtheta = rhs_CR3BP_3D(theta_rpo, X_rpo_post(:), nlp.mu_em);

X_rpo_post = X_rpo_post(:).';
dXrpo_dt = zeros(6, 1); % RPO boundary model depends on theta only in this strict path.

% Map regularized impulses to delta-v vectors and Jacobians
% Allocate per-node DV vectors and d(DV)/dU blocks.
DV_nodes = zeros(N, 3);
dDVdU_nodes = zeros(3, 3, N);

% Build each node DV mapping using case-2a or case-2b formula.
for j = 1:N
    u = U_nodes(j, 1);
    w = U_nodes(j, 2);
    s = U_nodes(j, 3);

    is_case2b = any(nlp.I_case2b(:) == j);

    if ~is_case2b
        DV_nodes(j, :) = [u^2 - w^2 - s^2, 2 * u * w, 2 * u * s];
        dDVdU_nodes(:, :, j) = [2 * u, -2 * w, -2 * s; ...
                                2 * w,  2 * u,       0; ...
                                2 * s,      0,   2 * u];
    else
        DV_nodes(j, :) = [2 * u * w, u^2 - w^2 - s^2, 2 * u * s];
        dDVdU_nodes(:, :, j) = [2 * w,  2 * u,       0; ...
                                2 * u, -2 * w, -2 * s; ...
                                2 * s,      0,   2 * u];
    end
end

% Reconstruct full node states
% Allocate all node states and fill interior direct decision states.
X_nodes = zeros(N, 6);
X_nodes(2:N-1, :) = X_nodes_interior(2:N-1, :);

% Build first node as SEL2 pre-state plus node-1 impulse on velocity.
X1_post = X_sel2_pre(:);
X1_post(4:6) = X1_post(4:6) + DV_nodes(1, :).';
X_nodes(1, :) = X1_post(:).';

% Build last node as phase-defined RPO post-state.
X_nodes(N, :) = X_rpo_post(:).';

% Constraint sizing and allocation
% Count inequality rows: dt positivity, always-on max TOF, and peri rows.
nc_time = N - 1;
nc_tof = 1;
nc_peri_ineq = 2 + 2 * double(use_earth_peri);
nc = nc_time + nc_tof + nc_peri_ineq;

% Count equality rows: defects plus Moon peri (and optional Earth peri).
neq_defects = 6 * (N - 1);
neq_peri = 1 + double(use_earth_peri);
neq = neq_defects + neq_peri;

% Allocate constraint vectors.
c = zeros(nc, 1);
ceq = zeros(neq, 1);

% Allocate sparse Jacobians when requested.
if nargout > 2
    Jc = spalloc(nc, nvars, max(64, 8 * nc));
    Jceq = spalloc(neq, nvars, max(256, 24 * neq));
else
    Jc = [];
    Jceq = [];
end

% Time inequalities
% Enforce strictly positive segment durations.
for i = 1:(N - 1)
    c(i) = nlp.min_dt_nd - dt_nodes(i);

    if nargout > 2
        Jc(i, nlp.map.dt_nodes(i)) = -1;
    end
end

% Enforce always-on transfer-time upper bound.
row_tof = nc_time + 1;
c(row_tof) = transfer_time_nd - nlp.max_tof_nd;

if nargout > 2
    Jc(row_tof, nlp.map.dt_nodes(:)) = 1;
end

% Multiple-shooting defects
% Build segment defects and Jacobian blocks for each arc i -> i+1.
for i = 1:(N - 1)
    % Read segment endpoint states and times.
    Xi = X_nodes(i, :).';
    Xnext = X_nodes(i + 1, :).';
    ti = t_nodes(i);
    tip1 = t_nodes(i + 1);

    % Propagate one BCR4BP segment with STM enabled.
    [Xprop, PHI] = propagate(Xi, ti, tip1, nlp.mu_em, '3D', 'BCR4BP', 'true', ...
        [], [], nlp.theta_sun0_rad, nlp.abs_tol, nlp.rel_tol);

    % Read arrival impulse on node i+1.
    DV_next = DV_nodes(i + 1, :).';

    % Define row indices for this segment.
    row0 = 6 * (i - 1);
    rows_pos = row0 + (1:3);
    rows_vel = row0 + (4:6);

    % Enforce position continuity defects.
    ceq(rows_pos) = Xprop(1:3) - Xnext(1:3);

    % Enforce velocity continuity defects with arrival impulse.
    ceq(rows_vel) = Xnext(4:6) - Xprop(4:6) - DV_next;

    if nargout > 2
        % Evaluate segment endpoint dynamics for time sensitivities.
        fi = rhs_BCR4BP_3D(ti, Xi, nlp.mu_em, nlp.theta_sun0_rad);
        ff = rhs_BCR4BP_3D(tip1, Xprop, nlp.mu_em, nlp.theta_sun0_rad);

        % Build propagated-state sensitivity wrt segment start and end times.
        dXprop_dti = -PHI * fi;
        dXprop_dtf = ff;

        % Add start-node state sensitivities.
        if i == 1
            % Map first-node state sensitivity wrt theta_sel2 and U1.
            dXi_dtheta_sel2 = dXsel2_dtheta(:);
            dXi_dt1 = dXsel2_dt(:);
            dXi_dU1 = [zeros(3, 3); dDVdU_nodes(:, :, 1)];

            % Chain theta_sel2 through segment propagation.
            Jceq(rows_pos, nlp.map.theta_sel2) = Jceq(rows_pos, nlp.map.theta_sel2) + PHI(1:3, :) * dXi_dtheta_sel2;
            Jceq(rows_vel, nlp.map.theta_sel2) = Jceq(rows_vel, nlp.map.theta_sel2) - PHI(4:6, :) * dXi_dtheta_sel2;

            % Chain U1 through segment propagation.
            idxU1 = nlp.map.node(1).U;
            Jceq(rows_pos, idxU1) = Jceq(rows_pos, idxU1) + PHI(1:3, :) * dXi_dU1;
            Jceq(rows_vel, idxU1) = Jceq(rows_vel, idxU1) - PHI(4:6, :) * dXi_dU1;

            % Add explicit t1 dependence from boundary-state time sensitivity.
            dXprop_dti = dXprop_dti + PHI * dXi_dt1;
        else
            % Map interior start-node state directly to decision X variables.
            idxXi = nlp.map.node(i).X;
            Jceq(rows_pos, idxXi) = Jceq(rows_pos, idxXi) + PHI(1:3, 1:6);
            Jceq(rows_vel, idxXi) = Jceq(rows_vel, idxXi) - PHI(4:6, 1:6);
        end

        % Build defect derivatives wrt endpoint times.
        dceq_dti_pos = dXprop_dti(1:3);
        dceq_dti_vel = -dXprop_dti(4:6);

        if i + 1 == N
            dceq_dtip1_pos = dXprop_dtf(1:3) - dXrpo_dt(1:3);
            dceq_dtip1_vel = -dXprop_dtf(4:6) + dXrpo_dt(4:6);
        else
            dceq_dtip1_pos = dXprop_dtf(1:3);
            dceq_dtip1_vel = -dXprop_dtf(4:6);
        end

        % Map endpoint-time terms to theta_sel2 through start-wait coupling.
        Jceq(rows_pos, nlp.map.theta_sel2) = Jceq(rows_pos, nlp.map.theta_sel2) + ...
            dt_start_em_dtheta * (dceq_dti_pos + dceq_dtip1_pos);
        Jceq(rows_vel, nlp.map.theta_sel2) = Jceq(rows_vel, nlp.map.theta_sel2) + ...
            dt_start_em_dtheta * (dceq_dti_vel + dceq_dtip1_vel);

        % Map d(t_i) and d(t_{i+1}) terms to cumulative dt decision variables.
        if i > 1
            idx_dt_before = nlp.map.dt_nodes(1:(i - 1));
            Jceq(rows_pos, idx_dt_before) = Jceq(rows_pos, idx_dt_before) + dceq_dti_pos * ones(1, numel(idx_dt_before));
            Jceq(rows_vel, idx_dt_before) = Jceq(rows_vel, idx_dt_before) + dceq_dti_vel * ones(1, numel(idx_dt_before));
        end

        idx_dt_up_to = nlp.map.dt_nodes(1:i);
        Jceq(rows_pos, idx_dt_up_to) = Jceq(rows_pos, idx_dt_up_to) + dceq_dtip1_pos * ones(1, numel(idx_dt_up_to));
        Jceq(rows_vel, idx_dt_up_to) = Jceq(rows_vel, idx_dt_up_to) + dceq_dtip1_vel * ones(1, numel(idx_dt_up_to));

        % Add arrival-node state sensitivities.
        if i + 1 < N
            idxXnext = nlp.map.node(i + 1).X;
            Jceq(rows_pos, idxXnext(1:3)) = Jceq(rows_pos, idxXnext(1:3)) - eye(3);
            Jceq(rows_vel, idxXnext(4:6)) = Jceq(rows_vel, idxXnext(4:6)) + eye(3);
        else
            Jceq(rows_pos, nlp.map.theta_rpo) = Jceq(rows_pos, nlp.map.theta_rpo) - dXrpo_dtheta(1:3);
            Jceq(rows_vel, nlp.map.theta_rpo) = Jceq(rows_vel, nlp.map.theta_rpo) + dXrpo_dtheta(4:6);
        end

        % Add arrival-impulse sensitivities in velocity rows.
        idxU_next = nlp.map.node(i + 1).U;
        Jceq(rows_vel, idxU_next) = Jceq(rows_vel, idxU_next) - dDVdU_nodes(:, :, i + 1);
    end
end

% Moon periapsis constraints
idx_peri = round(nlp.idx_peri_node);
X_peri = X_nodes(idx_peri, :).';
t_peri = t_nodes(idx_peri);

% Evaluate Moon-relative range, radial speed, and radial acceleration.
r_vec_m = [X_peri(1) - 1 + nlp.mu_em; X_peri(2); X_peri(3)];
v_vec_m = X_peri(4:6);
r2 = norm(r_vec_m);
f_peri_m = rhs_BCR4BP_3D(t_peri, X_peri, nlp.mu_em, nlp.theta_sun0_rad);
a_vec_m = f_peri_m(4:6);
rv_m = dot(r_vec_m, v_vec_m);
r2dot = rv_m / r2;
r2ddot = (dot(r_vec_m, a_vec_m) + dot(v_vec_m, v_vec_m)) / r2 - (rv_m^2) / (r2^3);

% Build Moon inequality and equality rows.
row_rmin_moon = nc_time + nc_tof + 1;
row_rddot_moon = nc_time + nc_tof + 2;
row_peri_eq_moon = neq_defects + 1;

c(row_rmin_moon) = nlp.rmin_moon_nd^2 - r2^2;
c(row_rddot_moon) = -r2ddot;
ceq(row_peri_eq_moon) = r2dot;

if nargout > 2
    % Build analytic gradient for Moon radius-floor inequality.
    grad_rmin_X = zeros(1, 6);
    grad_rmin_X(1) = -2 * (X_peri(1) - 1 + nlp.mu_em);
    grad_rmin_X(2) = -2 * X_peri(2);
    grad_rmin_X(3) = -2 * X_peri(3);

    % Build finite-difference gradients for Moon radial rate and acceleration.
    grad_r2dot_X = zeros(1, 6);
    grad_r2ddot_X = zeros(1, 6);

    for k = 1:6
        hX = sqrt(eps) * max(1, abs(X_peri(k)));

        Xp = X_peri;
        Xm = X_peri;
        Xp(k) = Xp(k) + hX;
        Xm(k) = Xm(k) - hX;

        % Evaluate Moon kinematics at positive perturbation.
        r_vec_p = [Xp(1) - 1 + nlp.mu_em; Xp(2); Xp(3)];
        v_vec_p = Xp(4:6);
        r_p = norm(r_vec_p);
        f_p = rhs_BCR4BP_3D(t_peri, Xp, nlp.mu_em, nlp.theta_sun0_rad);
        a_vec_p = f_p(4:6);
        rv_p = dot(r_vec_p, v_vec_p);
        rdot_p = rv_p / r_p;
        rddot_p = (dot(r_vec_p, a_vec_p) + dot(v_vec_p, v_vec_p)) / r_p - (rv_p^2) / (r_p^3);

        % Evaluate Moon kinematics at negative perturbation.
        r_vec_mn = [Xm(1) - 1 + nlp.mu_em; Xm(2); Xm(3)];
        v_vec_mn = Xm(4:6);
        r_mn = norm(r_vec_mn);
        f_mn = rhs_BCR4BP_3D(t_peri, Xm, nlp.mu_em, nlp.theta_sun0_rad);
        a_vec_mn = f_mn(4:6);
        rv_mn = dot(r_vec_mn, v_vec_mn);
        rdot_mn = rv_mn / r_mn;
        rddot_mn = (dot(r_vec_mn, a_vec_mn) + dot(v_vec_mn, v_vec_mn)) / r_mn - (rv_mn^2) / (r_mn^3);

        grad_r2dot_X(k) = (rdot_p - rdot_mn) / (2 * hX);
        grad_r2ddot_X(k) = (rddot_p - rddot_mn) / (2 * hX);
    end

    % Build finite-difference time gradients for Moon radial rate and acceleration.
    ht_peri = sqrt(eps) * max(1, abs(t_peri));

    f_tp = rhs_BCR4BP_3D(t_peri + ht_peri, X_peri, nlp.mu_em, nlp.theta_sun0_rad);
    a_tp = f_tp(4:6);
    rdot_tp = rv_m / r2;
    rddot_tp = (dot(r_vec_m, a_tp) + dot(v_vec_m, v_vec_m)) / r2 - (rv_m^2) / (r2^3);

    f_tm = rhs_BCR4BP_3D(t_peri - ht_peri, X_peri, nlp.mu_em, nlp.theta_sun0_rad);
    a_tm = f_tm(4:6);
    rdot_tm = rv_m / r2;
    rddot_tm = (dot(r_vec_m, a_tm) + dot(v_vec_m, v_vec_m)) / r2 - (rv_m^2) / (r2^3);

    grad_r2dot_t = (rdot_tp - rdot_tm) / (2 * ht_peri);
    grad_r2ddot_t = (rddot_tp - rddot_tm) / (2 * ht_peri);

    % Route Moon radius-floor gradient from node-state space to decision space.
    if idx_peri > 1 && idx_peri < N
        idxX_peri = nlp.map.node(idx_peri).X;
        Jc(row_rmin_moon, idxX_peri) = Jc(row_rmin_moon, idxX_peri) + grad_rmin_X;
    elseif idx_peri == 1
        dX1_dtheta_sel2 = dXsel2_dtheta(:);
        dX1_dU1 = [zeros(3, 3); dDVdU_nodes(:, :, 1)];
        idxU1 = nlp.map.node(1).U;

        Jc(row_rmin_moon, nlp.map.theta_sel2) = Jc(row_rmin_moon, nlp.map.theta_sel2) + grad_rmin_X * dX1_dtheta_sel2;
        Jc(row_rmin_moon, idxU1) = Jc(row_rmin_moon, idxU1) + grad_rmin_X * dX1_dU1;
    else
        Jc(row_rmin_moon, nlp.map.theta_rpo) = Jc(row_rmin_moon, nlp.map.theta_rpo) + grad_rmin_X * dXrpo_dtheta(:);
    end

    % Route Moon peri-equality gradient from node-state space to decision space.
    if idx_peri > 1 && idx_peri < N
        idxX_peri = nlp.map.node(idx_peri).X;
        Jceq(row_peri_eq_moon, idxX_peri) = Jceq(row_peri_eq_moon, idxX_peri) + grad_r2dot_X;
    elseif idx_peri == 1
        dX1_dtheta_sel2 = dXsel2_dtheta(:);
        dX1_dU1 = [zeros(3, 3); dDVdU_nodes(:, :, 1)];
        idxU1 = nlp.map.node(1).U;

        Jceq(row_peri_eq_moon, nlp.map.theta_sel2) = Jceq(row_peri_eq_moon, nlp.map.theta_sel2) + grad_r2dot_X * dX1_dtheta_sel2;
        Jceq(row_peri_eq_moon, idxU1) = Jceq(row_peri_eq_moon, idxU1) + grad_r2dot_X * dX1_dU1;
    else
        Jceq(row_peri_eq_moon, nlp.map.theta_rpo) = Jceq(row_peri_eq_moon, nlp.map.theta_rpo) + grad_r2dot_X * dXrpo_dtheta(:);
    end

    % Route Moon radial-acceleration inequality gradient from state to decision space.
    grad_minus_r2ddot_X = -grad_r2ddot_X;
    if idx_peri > 1 && idx_peri < N
        idxX_peri = nlp.map.node(idx_peri).X;
        Jc(row_rddot_moon, idxX_peri) = Jc(row_rddot_moon, idxX_peri) + grad_minus_r2ddot_X;
    elseif idx_peri == 1
        dX1_dtheta_sel2 = dXsel2_dtheta(:);
        dX1_dU1 = [zeros(3, 3); dDVdU_nodes(:, :, 1)];
        idxU1 = nlp.map.node(1).U;

        Jc(row_rddot_moon, nlp.map.theta_sel2) = Jc(row_rddot_moon, nlp.map.theta_sel2) + grad_minus_r2ddot_X * dX1_dtheta_sel2;
        Jc(row_rddot_moon, idxU1) = Jc(row_rddot_moon, idxU1) + grad_minus_r2ddot_X * dX1_dU1;
    else
        Jc(row_rddot_moon, nlp.map.theta_rpo) = Jc(row_rddot_moon, nlp.map.theta_rpo) + grad_minus_r2ddot_X * dXrpo_dtheta(:);
    end

    % Route Moon peri-equality direct time gradient to theta_sel2 and dt nodes.
    Jceq(row_peri_eq_moon, nlp.map.theta_sel2) = Jceq(row_peri_eq_moon, nlp.map.theta_sel2) + grad_r2dot_t * dt_start_em_dtheta;
    if idx_peri > 1
        idx_dt_peri = nlp.map.dt_nodes(1:(idx_peri - 1));
        Jceq(row_peri_eq_moon, idx_dt_peri) = Jceq(row_peri_eq_moon, idx_dt_peri) + grad_r2dot_t;
    end

    % Route Moon radial-acceleration direct time gradient to theta_sel2 and dt nodes.
    Jc(row_rddot_moon, nlp.map.theta_sel2) = Jc(row_rddot_moon, nlp.map.theta_sel2) - grad_r2ddot_t * dt_start_em_dtheta;
    if idx_peri > 1
        idx_dt_peri = nlp.map.dt_nodes(1:(idx_peri - 1));
        Jc(row_rddot_moon, idx_dt_peri) = Jc(row_rddot_moon, idx_dt_peri) - grad_r2ddot_t;
    end
end

% Optional Earth periapsis constraints
if use_earth_peri
    idx_earth_peri = round(nlp.idx_earth_peri_node);
    X_earth_peri = X_nodes(idx_earth_peri, :).';
    t_earth_peri = t_nodes(idx_earth_peri);

    % Evaluate Earth-relative range, radial speed, and radial acceleration.
    r_vec_e = [X_earth_peri(1) + nlp.mu_em; X_earth_peri(2); X_earth_peri(3)];
    v_vec_e = X_earth_peri(4:6);
    r1 = norm(r_vec_e);
    f_peri_e = rhs_BCR4BP_3D(t_earth_peri, X_earth_peri, nlp.mu_em, nlp.theta_sun0_rad);
    a_vec_e = f_peri_e(4:6);
    rv_e = dot(r_vec_e, v_vec_e);
    r1dot = rv_e / r1;
    r1ddot = (dot(r_vec_e, a_vec_e) + dot(v_vec_e, v_vec_e)) / r1 - (rv_e^2) / (r1^3);

    % Build Earth inequality and equality rows.
    row_rmin_earth = nc_time + nc_tof + 3;
    row_rddot_earth = nc_time + nc_tof + 4;
    row_peri_eq_earth = neq_defects + 2;

    c(row_rmin_earth) = nlp.rmin_earth_nd^2 - r1^2;
    c(row_rddot_earth) = -r1ddot;
    ceq(row_peri_eq_earth) = r1dot;

    if nargout > 2
        % Build analytic gradient for Earth radius-floor inequality.
        grad_rmin_earth_X = zeros(1, 6);
        grad_rmin_earth_X(1) = -2 * (X_earth_peri(1) + nlp.mu_em);
        grad_rmin_earth_X(2) = -2 * X_earth_peri(2);
        grad_rmin_earth_X(3) = -2 * X_earth_peri(3);

        % Build finite-difference gradients for Earth radial rate and acceleration.
        grad_r1dot_X = zeros(1, 6);
        grad_r1ddot_X = zeros(1, 6);

        for k = 1:6
            hX = sqrt(eps) * max(1, abs(X_earth_peri(k)));

            Xp = X_earth_peri;
            Xm = X_earth_peri;
            Xp(k) = Xp(k) + hX;
            Xm(k) = Xm(k) - hX;

            % Evaluate Earth kinematics at positive perturbation.
            r_vec_p = [Xp(1) + nlp.mu_em; Xp(2); Xp(3)];
            v_vec_p = Xp(4:6);
            r_p = norm(r_vec_p);
            f_p = rhs_BCR4BP_3D(t_earth_peri, Xp, nlp.mu_em, nlp.theta_sun0_rad);
            a_vec_p = f_p(4:6);
            rv_p = dot(r_vec_p, v_vec_p);
            rdot_p = rv_p / r_p;
            rddot_p = (dot(r_vec_p, a_vec_p) + dot(v_vec_p, v_vec_p)) / r_p - (rv_p^2) / (r_p^3);

            % Evaluate Earth kinematics at negative perturbation.
            r_vec_mn = [Xm(1) + nlp.mu_em; Xm(2); Xm(3)];
            v_vec_mn = Xm(4:6);
            r_mn = norm(r_vec_mn);
            f_mn = rhs_BCR4BP_3D(t_earth_peri, Xm, nlp.mu_em, nlp.theta_sun0_rad);
            a_vec_mn = f_mn(4:6);
            rv_mn = dot(r_vec_mn, v_vec_mn);
            rdot_mn = rv_mn / r_mn;
            rddot_mn = (dot(r_vec_mn, a_vec_mn) + dot(v_vec_mn, v_vec_mn)) / r_mn - (rv_mn^2) / (r_mn^3);

            grad_r1dot_X(k) = (rdot_p - rdot_mn) / (2 * hX);
            grad_r1ddot_X(k) = (rddot_p - rddot_mn) / (2 * hX);
        end

        % Build finite-difference time gradients for Earth radial rate and acceleration.
        ht_peri = sqrt(eps) * max(1, abs(t_earth_peri));

        f_tp = rhs_BCR4BP_3D(t_earth_peri + ht_peri, X_earth_peri, nlp.mu_em, nlp.theta_sun0_rad);
        a_tp = f_tp(4:6);
        rdot_tp = rv_e / r1;
        rddot_tp = (dot(r_vec_e, a_tp) + dot(v_vec_e, v_vec_e)) / r1 - (rv_e^2) / (r1^3);

        f_tm = rhs_BCR4BP_3D(t_earth_peri - ht_peri, X_earth_peri, nlp.mu_em, nlp.theta_sun0_rad);
        a_tm = f_tm(4:6);
        rdot_tm = rv_e / r1;
        rddot_tm = (dot(r_vec_e, a_tm) + dot(v_vec_e, v_vec_e)) / r1 - (rv_e^2) / (r1^3);

        grad_r1dot_t = (rdot_tp - rdot_tm) / (2 * ht_peri);
        grad_r1ddot_t = (rddot_tp - rddot_tm) / (2 * ht_peri);

        % Route Earth radius-floor gradient from node-state space to decision space.
        if idx_earth_peri > 1 && idx_earth_peri < N
            idxX_earth_peri = nlp.map.node(idx_earth_peri).X;
            Jc(row_rmin_earth, idxX_earth_peri) = Jc(row_rmin_earth, idxX_earth_peri) + grad_rmin_earth_X;
        elseif idx_earth_peri == 1
            dX1_dtheta_sel2 = dXsel2_dtheta(:);
            dX1_dU1 = [zeros(3, 3); dDVdU_nodes(:, :, 1)];
            idxU1 = nlp.map.node(1).U;

            Jc(row_rmin_earth, nlp.map.theta_sel2) = Jc(row_rmin_earth, nlp.map.theta_sel2) + grad_rmin_earth_X * dX1_dtheta_sel2;
            Jc(row_rmin_earth, idxU1) = Jc(row_rmin_earth, idxU1) + grad_rmin_earth_X * dX1_dU1;
        else
            Jc(row_rmin_earth, nlp.map.theta_rpo) = Jc(row_rmin_earth, nlp.map.theta_rpo) + grad_rmin_earth_X * dXrpo_dtheta(:);
        end

        % Route Earth peri-equality gradient from node-state space to decision space.
        if idx_earth_peri > 1 && idx_earth_peri < N
            idxX_earth_peri = nlp.map.node(idx_earth_peri).X;
            Jceq(row_peri_eq_earth, idxX_earth_peri) = Jceq(row_peri_eq_earth, idxX_earth_peri) + grad_r1dot_X;
        elseif idx_earth_peri == 1
            dX1_dtheta_sel2 = dXsel2_dtheta(:);
            dX1_dU1 = [zeros(3, 3); dDVdU_nodes(:, :, 1)];
            idxU1 = nlp.map.node(1).U;

            Jceq(row_peri_eq_earth, nlp.map.theta_sel2) = Jceq(row_peri_eq_earth, nlp.map.theta_sel2) + grad_r1dot_X * dX1_dtheta_sel2;
            Jceq(row_peri_eq_earth, idxU1) = Jceq(row_peri_eq_earth, idxU1) + grad_r1dot_X * dX1_dU1;
        else
            Jceq(row_peri_eq_earth, nlp.map.theta_rpo) = Jceq(row_peri_eq_earth, nlp.map.theta_rpo) + grad_r1dot_X * dXrpo_dtheta(:);
        end

        % Route Earth radial-acceleration inequality gradient from state to decision space.
        grad_minus_r1ddot_X = -grad_r1ddot_X;
        if idx_earth_peri > 1 && idx_earth_peri < N
            idxX_earth_peri = nlp.map.node(idx_earth_peri).X;
            Jc(row_rddot_earth, idxX_earth_peri) = Jc(row_rddot_earth, idxX_earth_peri) + grad_minus_r1ddot_X;
        elseif idx_earth_peri == 1
            dX1_dtheta_sel2 = dXsel2_dtheta(:);
            dX1_dU1 = [zeros(3, 3); dDVdU_nodes(:, :, 1)];
            idxU1 = nlp.map.node(1).U;

            Jc(row_rddot_earth, nlp.map.theta_sel2) = Jc(row_rddot_earth, nlp.map.theta_sel2) + grad_minus_r1ddot_X * dX1_dtheta_sel2;
            Jc(row_rddot_earth, idxU1) = Jc(row_rddot_earth, idxU1) + grad_minus_r1ddot_X * dX1_dU1;
        else
            Jc(row_rddot_earth, nlp.map.theta_rpo) = Jc(row_rddot_earth, nlp.map.theta_rpo) + grad_minus_r1ddot_X * dXrpo_dtheta(:);
        end

        % Route Earth peri-equality direct time gradient to theta_sel2 and dt nodes.
        Jceq(row_peri_eq_earth, nlp.map.theta_sel2) = Jceq(row_peri_eq_earth, nlp.map.theta_sel2) + grad_r1dot_t * dt_start_em_dtheta;
        if idx_earth_peri > 1
            idx_dt_earth_peri = nlp.map.dt_nodes(1:(idx_earth_peri - 1));
            Jceq(row_peri_eq_earth, idx_dt_earth_peri) = Jceq(row_peri_eq_earth, idx_dt_earth_peri) + grad_r1dot_t;
        end

        % Route Earth radial-acceleration direct time gradient to theta_sel2 and dt nodes.
        Jc(row_rddot_earth, nlp.map.theta_sel2) = Jc(row_rddot_earth, nlp.map.theta_sel2) - grad_r1ddot_t * dt_start_em_dtheta;
        if idx_earth_peri > 1
            idx_dt_earth_peri = nlp.map.dt_nodes(1:(idx_earth_peri - 1));
            Jc(row_rddot_earth, idx_dt_earth_peri) = Jc(row_rddot_earth, idx_dt_earth_peri) - grad_r1ddot_t;
        end
    end
end

% Return Jacobians in fmincon orientation
if nargout > 2
    Jc = Jc.';
    Jceq = Jceq.';
end
end
