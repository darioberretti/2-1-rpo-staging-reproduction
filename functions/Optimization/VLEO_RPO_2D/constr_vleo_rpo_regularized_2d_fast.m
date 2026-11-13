function [c, ceq, Jc, Jceq] = constr_vleo_rpo_regularized_2d_fast(y, nlp)
%CONSTR_VLEO_RPO_REGULARIZED_2D_FAST Evaluate the planar BCR4BP shooting constraints and sparse Jacobians.
%
% PROTOTYPE:
%   [c, ceq, Jc, Jceq] = constr_vleo_rpo_regularized_2d_fast(y, nlp)
%
% DESCRIPTION:
%   Build the first post-departure and last pre-insertion states from the
%   circular parking orbit and free-phase CR3BP RPO target. Enforce segment
%   continuity, node-wise Moon clearance, and the selected impulse/altitude/
%   TOF policies. Integrate each segment with STM and initial-Sun-angle
%   sensitivity when derivatives are requested.
%
% INPUT:
%   y    [nvars x 1] Mapped vector from build_vleo_rpo_ms_initial_guess_2d.
%   nlp  setup.nlp with node map, normalized model/scales, departure/arrival
%          boundaries, impulse branches, and timing/clearance policies.
%          States use DU/VU, segment durations and RPO phase TU, departure/
%          Sun phases rad, optional departure altitude km.
%
% OUTPUT:
%   c, ceq    Columns feasible when c <= 0 and ceq = 0. Moon-distance
%               inequalities apply at nodes, not continuously along each arc.
%   Jc, Jceq  Optional sparse [nvars x number of constraints] Jacobians in
%               fmincon convention. Residual rows retain their native units.

if nargout > 2
    [c, ceq, Jc, Jceq] = local_evaluate_constraints_fast(y, nlp, true);

    % Transpose the Jacobians into the fmincon convention.
    Jc = Jc.';
    Jceq = Jceq.';
else
    % Evaluate only the constraint values when no Jacobian is requested.
    [c, ceq] = local_evaluate_constraints_fast(y, nlp, false);
end

end

function [c, ceq, Jc, Jceq] = local_evaluate_constraints_fast(y, nlp, compute_jacobian)
%LOCAL_EVALUATE_CONSTRAINTS_FAST Assemble endpoint states, segment defects, and active policy rows.
%
% PROTOTYPE:
%   [c, ceq, Jc, Jceq] = local_evaluate_constraints_fast(y, nlp, compute_jacobian)
%
% DESCRIPTION:
%   Decode the mapped state/boundary variables, propagate each segment, and
%   assemble continuity plus impulse/altitude/TOF/Moon-distance residuals.
%   Keep internal Jacobians in row-per-constraint orientation; the public
%   wrapper transposes them for fmincon.
%
% INPUT:
%   y                 Mapped decision-vector column.
%   nlp               Planar problem schema documented by the public wrapper.
%   compute_jacobian  Logical enabling state/time/boundary sensitivities.
%
% OUTPUT:
%   c, ceq            Inequality/equality residual columns.
%   Jc, Jceq          Sparse [number of constraints x numel(y)] matrices,
%                       or [] when sensitivities are disabled.

dec = decode_vleo_rpo_decision_vector_2d(y, nlp);

% Read the node count and the node times.
N = nlp.N;
t_nodes = dec.t_nodes(:);
n_vars = numel(y);

% Build the departure-circle state and its derivatives.
[x_departure_pre, dxdep_dphase, dxdep_dh] = circular_parking_state_2d(dec.phase_departure, dec.departure_altitude_km, nlp);

% Map the departure regularized variables into a planar delta-v.
[dv_departure, dDV1dU] = planar_oshima_map_2d(dec.U_nodes(1, :).', any(nlp.I_case2b == 1));

% Build the first post-impulse node state.
x_first = x_departure_pre;
x_first(3:4) = x_first(3:4) + dv_departure;

% Build the free-phase nominal 3BP RPO target state.
[x_rpo_target, dxrpo_dphase] = local_rpo_target_state_fast(dec.phase_rpo, nlp);

% Map the arrival regularized variables into a planar delta-v.
[dv_arrival, dDVNdU] = planar_oshima_map_2d(dec.U_nodes(end, :).', any(nlp.I_case2b == N));

% Build the last pre-impulse node state.
x_last = x_rpo_target;
x_last(3:4) = x_last(3:4) - dv_arrival;

% Assemble the full node-state array.
X_nodes = zeros(N, 4);
X_nodes(1, :) = x_first(:).';
X_nodes(end, :) = x_last(:).';
if N > 2
    X_nodes(2:N-1, :) = dec.X_nodes_interior(2:N-1, :);
end

% Count the inequality rows.
n_time_ineq = 0;
n_tof_ineq = double(strcmpi(nlp.tof_mode, 'free_with_max_bound'));
n_moon_ineq = N;
n_ineq = n_time_ineq + n_tof_ineq + n_moon_ineq;

% Count the equality rows.
n_defect_eq = 4 * (N - 1);
n_zero_u_eq = 0;
if strcmpi(nlp.impulse_mode, 'endpoints_only')
    n_zero_u_eq = 2 * max(0, N - 2);
end
n_altitude_eq = double(strcmpi(nlp.departure_constraint_mode, 'fixed_altitude_seed_warmstart'));
n_tof_eq = double(strcmpi(nlp.tof_mode, 'fixed_to_seed'));
n_eq = n_defect_eq + n_zero_u_eq + n_altitude_eq + n_tof_eq;

% Allocate the constraint vectors.
c = zeros(n_ineq, 1);
ceq = zeros(n_eq, 1);

% Allocate the sparse Jacobians only when requested.
if compute_jacobian
    Jc = spalloc(n_ineq, n_vars, max(16, 4 * N));
    Jceq = spalloc(n_eq, n_vars, max(64, 28 * (N - 1) + 8 * N));
else
    Jc = [];
    Jceq = [];
end

% Build the multiple-shooting defects one segment at a time.
for i = 1:(N - 1)
    % Read the segment endpoints and times.
    Xi = X_nodes(i, :).';
    Xnext = X_nodes(i + 1, :).';
    ti = t_nodes(i);
    tf = t_nodes(i + 1);

    % Propagate the current segment with or without sensitivities.
    if compute_jacobian
        [Xprop, PHI, dXprop_dtheta] = propagate_bcr4bp_2d_stm_theta_fast( ...
            Xi, ti, tf, nlp.mu, dec.theta_sun0, nlp.abs_tol, nlp.rel_tol, true);
    else
        Xprop = propagate_bcr4bp_2d_stm_theta_fast( ...
            Xi, ti, tf, nlp.mu, dec.theta_sun0, nlp.abs_tol, nlp.rel_tol, false);
        PHI = [];
        dXprop_dtheta = [];
    end

    % Store the segment defect.
    row_defect = 4 * (i - 1) + (1:4);
    ceq(row_defect) = Xprop - Xnext;

    % Differentiate the defect analytically when requested.
    if compute_jacobian
        % Evaluate the segment endpoint dynamics for time derivatives.
        fi = rhs_BCR4BP_2D(ti, Xi, nlp.mu, dec.theta_sun0);
        ff = rhs_BCR4BP_2D(tf, Xprop, nlp.mu, dec.theta_sun0);

        % Build the propagated-state derivatives with respect to endpoint times.
        dXprop_dti = -PHI * fi;
        dXprop_dtf = ff;

        % Route the start-node derivatives.
        if i == 1
            % Route the departure-phase derivative through the first node.
            Jceq(row_defect, nlp.map.phase_departure) = Jceq(row_defect, nlp.map.phase_departure) + PHI * dxdep_dphase;

            % Route the departure-altitude derivative only when it exists.
            if ~isempty(nlp.map.departure_altitude_km)
                Jceq(row_defect, nlp.map.departure_altitude_km) = Jceq(row_defect, nlp.map.departure_altitude_km) + PHI * dxdep_dh;
            end

            % Route the first-node impulse derivative through the first node.
            Jceq(row_defect, nlp.map.node(1).U) = Jceq(row_defect, nlp.map.node(1).U) + PHI(:, 3:4) * dDV1dU;
        else
            % Route the interior start-node state derivative directly.
            Jceq(row_defect, nlp.map.node(i).X) = Jceq(row_defect, nlp.map.node(i).X) + PHI;
        end

        % Route the end-node derivatives.
        if i + 1 < N
            % Route the interior end-node state derivative directly.
            Jceq(row_defect, nlp.map.node(i + 1).X) = Jceq(row_defect, nlp.map.node(i + 1).X) - speye(4);
        else
            % Route the arrival-phase derivative through the nominal RPO target.
            Jceq(row_defect, nlp.map.phase_rpo) = Jceq(row_defect, nlp.map.phase_rpo) - dxrpo_dphase;

            % Route the arrival impulse derivative through the last node.
            Jceq(row_defect(3:4), nlp.map.node(N).U) = Jceq(row_defect(3:4), nlp.map.node(N).U) + dDVNdU;
        end

        % Route the segment-duration derivatives.
        if i > 1
            idx_dt_before = nlp.map.dt_nodes(1:(i - 1));
            Jceq(row_defect, idx_dt_before) = Jceq(row_defect, idx_dt_before) + (dXprop_dti + dXprop_dtf) * ones(1, numel(idx_dt_before));
        end
        Jceq(row_defect, nlp.map.dt_nodes(i)) = Jceq(row_defect, nlp.map.dt_nodes(i)) + dXprop_dtf;

        % Route the direct thetaSun0 sensitivity from the segment propagation.
        Jceq(row_defect, nlp.map.theta_sun0) = Jceq(row_defect, nlp.map.theta_sun0) + dXprop_dtheta;
    end
end

% Enforce zero interior impulses in endpoint-only mode.
if strcmpi(nlp.impulse_mode, 'endpoints_only') && N > 2
    % Start writing after the defect rows.
    row_zero = n_defect_eq;

    % Enforce one zero-impulse condition per interior node component.
    for j = 2:(N - 1)
        % Route the interior regularized variables directly.
        ceq(row_zero + 1:row_zero + 2) = dec.U_nodes(j, :).';

        % Differentiate the zero-impulse constraint directly.
        if compute_jacobian
            Jceq(row_zero + 1:row_zero + 2, nlp.map.node(j).U) = speye(2);
        end

        % Advance the equality-row counter.
        row_zero = row_zero + 2;
    end
end

% Start the scalar equality row cursor after the node-defect and zero-impulse rows.
row_scalar_eq = n_defect_eq + n_zero_u_eq;

% Enforce the requested fixed altitude after warm-starting from the seed altitude.
if strcmpi(nlp.departure_constraint_mode, 'fixed_altitude_seed_warmstart')
    % Advance to the fixed-altitude equality row.
    row_scalar_eq = row_scalar_eq + 1;

    % Store the altitude mismatch scaled by the Earth-Moon distance unit.
    ceq(row_scalar_eq) = (dec.departure_altitude_km - nlp.departure_altitude_target_km) / nlp.DU;

    % Differentiate the fixed-altitude equality directly.
    if compute_jacobian
        Jceq(row_scalar_eq, nlp.map.departure_altitude_km) = 1 / nlp.DU;
    end
end

% Enforce a fixed total time of flight when requested.
if strcmpi(nlp.tof_mode, 'fixed_to_seed')
    % Advance to the fixed-TOF equality row.
    row_scalar_eq = row_scalar_eq + 1;

    % Store the total-time equality after the other scalar equalities.
    ceq(row_scalar_eq) = sum(dec.dt_nodes) - nlp.seed_tof_tu;

    % Differentiate the total-time equality directly.
    if compute_jacobian
        Jceq(row_scalar_eq, nlp.map.dt_nodes) = 1;
    end
end

% Enforce the upper time-of-flight bound when requested.
if strcmpi(nlp.tof_mode, 'free_with_max_bound')
    % Store the total-time inequality at the start of the inequality vector.
    c(1) = sum(dec.dt_nodes) - nlp.max_tof_tu;

    % Differentiate the total-time inequality directly.
    if compute_jacobian
        Jc(1, nlp.map.dt_nodes) = 1;
    end
end

% Offset the Moon-distance rows when the TOF inequality exists.
moon_row_offset = n_tof_ineq;

% Enforce the Moon minimum-altitude condition on every node.
for j = 1:N
    % Read the planar Moon-relative position.
    r_moon = [X_nodes(j, 1) - 1 + nlp.mu; X_nodes(j, 2)];

    % Store the squared-distance inequality at node j.
    row_moon = moon_row_offset + j;
    c(row_moon) = nlp.rmin_moon_du^2 - dot(r_moon, r_moon);

    % Differentiate the node-wise Moon inequality when requested.
    if compute_jacobian
        % Build the planar position gradient.
        grad_pos = [-2 * r_moon(1), -2 * r_moon(2), 0, 0];

        % Route the first-node Moon gradient through the departure state.
        if j == 1
            Jc(row_moon, nlp.map.phase_departure) = Jc(row_moon, nlp.map.phase_departure) + grad_pos * dxdep_dphase;
            if ~isempty(nlp.map.departure_altitude_km)
                Jc(row_moon, nlp.map.departure_altitude_km) = Jc(row_moon, nlp.map.departure_altitude_km) + grad_pos * dxdep_dh;
            end

        % Route each interior-node Moon gradient directly.
        elseif j < N
            Jc(row_moon, nlp.map.node(j).X) = Jc(row_moon, nlp.map.node(j).X) + grad_pos;

        % Route the last-node Moon gradient through the RPO target.
        else
            Jc(row_moon, nlp.map.phase_rpo) = Jc(row_moon, nlp.map.phase_rpo) + grad_pos * dxrpo_dphase;
        end
    end
end

end

function [x_rpo_target, dxrpo_dphase] = local_rpo_target_state_fast(phase_rpo, nlp)
%LOCAL_RPO_TARGET_STATE_FAST Recover the nominal CR3BP target state and its time-phase derivative.
%
% PROTOTYPE:
%   [x_rpo_target, dxrpo_dphase] = local_rpo_target_state_fast(phase_rpo, nlp)
%
% DESCRIPTION:
%   Propagate the canonical planar RPO to the selected phase and evaluate
%   its CR3BP dynamics to obtain the derivative with respect to that phase.
%
% INPUT:
%   phase_rpo  Target orbit time offset [TU].
%   nlp        Canonical RPO state, mu, and abs_tol/rel_tol settings.
%
% OUTPUT:
%   x_rpo_target   [4 x 1] Rotating target state in DU/VU.
%   dxrpo_dphase   [4 x 1] State derivative with respect to phase [per TU].

x_rpo_target = propagate(nlp.x0_rpo_planar_3bp, 0, phase_rpo, nlp.mu, '2D', 'CR3BP', 'false', false, [], [], nlp.abs_tol, nlp.rel_tol);

% Differentiate the target state with respect to the target phase.
dxrpo_dphase = rhs_CR3BP_2D(phase_rpo, x_rpo_target, nlp.mu);

end
