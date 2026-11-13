function solution = postprocess_vleo_rpo_optimization_2d(y, setup, tag)
%POSTPROCESS_VLEO_RPO_OPTIMIZATION_2D Decode, propagate, and summarize one planar transfer vector.
%
% PROTOTYPE:
%   solution = postprocess_vleo_rpo_optimization_2d(y, setup, tag)
%
% DESCRIPTION:
%   Rebuild the phase-defined endpoint states and regularized burns, integrate
%   each BCR4BP segment, and report physical maneuver/arrival metrics. Evaluate
%   constraint residuals with the legacy constraint function. Minimum lunar
%   altitude is measured on the returned trajectory samples, not by a separate
%   continuous-time minimization.
%
% INPUT:
%   y      Mapped initial/optimized decision-vector column.
%   setup  Output from build_vleo_rpo_ms_initial_guess_2d, with nlp/constants.
%   tag    Optional descriptive label, default solution.
%
% OUTPUT:
%   solution  Struct containing y/dec, X_nodes [N x 4], t_traj [TU],
%               X_traj [Nt x 4] in rotating DU/VU, normalized impulse vectors
%               and magnitudes, physical total/boundary burns [m/s], TOF [days],
%               sampled Moon altitude [km], and equality/inequality residual
%               summaries. Stored inertial boundary speeds use km/s.

if nargin < 3
    tag = 'solution';
end

% Read the NLP description and decode the decision vector.
nlp = setup.nlp;
dec = decode_vleo_rpo_decision_vector_2d(y, nlp);

% Build the departure boundary state and its initial impulse.
[x_departure_pre, ~, ~] = circular_parking_state_2d(dec.phase_departure, dec.departure_altitude_km, nlp);
[dv_departure, ~] = planar_oshima_map_2d(dec.U_nodes(1, :).', any(nlp.I_case2b == 1));
x_first = x_departure_pre;
x_first(3:4) = x_first(3:4) + dv_departure;

% Build the arrival boundary state and its final impulse.
x_rpo_target = propagate(nlp.x0_rpo_planar_3bp, 0, dec.phase_rpo, nlp.mu, '2D', 'CR3BP', 'false', false, [], [], nlp.abs_tol, nlp.rel_tol);
[dv_arrival, ~] = planar_oshima_map_2d(dec.U_nodes(end, :).', any(nlp.I_case2b == nlp.N));
x_last = x_rpo_target;
x_last(3:4) = x_last(3:4) - dv_arrival;

% Convert the departure pre- and post-impulse states into inertial coordinates.
x_departure_pre_in = rot2in(x_departure_pre.', 0, nlp.mu);
x_departure_post_in = rot2in(x_first.', 0, nlp.mu);

% Convert the arrival pre-insertion and target states into inertial coordinates.
t_final = dec.t_nodes(end);
x_arrival_pre_in = rot2in(x_last.', t_final, nlp.mu);
x_arrival_target_in = rot2in(x_rpo_target.', t_final, nlp.mu);

% Build the Moon inertial state at the final node epoch.
x_moon_final_in = rot2in([1 - nlp.mu, 0, 0, 0], t_final, nlp.mu);

% Assemble the full node-state array.
X_nodes = zeros(nlp.N, 4);
X_nodes(1, :) = x_first(:).';
X_nodes(end, :) = x_last(:).';
if nlp.N > 2
    X_nodes(2:nlp.N-1, :) = dec.X_nodes_interior(2:nlp.N-1, :);
end

% Preallocate one cell entry per propagated segment.
t_segments = cell(nlp.N - 1, 1);
X_segments = cell(nlp.N - 1, 1);

% Propagate every segment to rebuild the full continuous trajectory.
for i = 1:(nlp.N - 1)
    % Propagate one segment on the optimized time grid.
    [~, ~, ~, x_seg, t_seg] = propagate(X_nodes(i, :).', dec.t_nodes(i), dec.t_nodes(i + 1), nlp.mu, '2D', 'BCR4BP', 'false', [], [], dec.theta_sun0, nlp.abs_tol, nlp.rel_tol);

    % Remove duplicate boundary samples after the first segment.
    if i > 1
        x_seg = x_seg(2:end, :);
        t_seg = t_seg(2:end, :);
    end

    % Store the current propagated segment.
    X_segments{i} = x_seg;
    t_segments{i} = t_seg;
end

% Concatenate the propagated segments once at the end.
X_traj = vertcat(X_segments{:});
t_traj = vertcat(t_segments{:});

% Convert every node impulse into its planar delta-v vector.
DV_nodes_rot = zeros(nlp.N, 2);
DV_mag_nodes_vu = zeros(nlp.N, 1);
for j = 1:nlp.N
    % Map the regularized variables at node j into planar delta-v.
    [dv_j, ~] = planar_oshima_map_2d(dec.U_nodes(j, :).', any(nlp.I_case2b == j));
    DV_nodes_rot(j, :) = dv_j(:).';
    DV_mag_nodes_vu(j) = norm(dv_j);
end

% Read the active node-wise delta-v magnitudes.
DV_mag_active_vu = DV_mag_nodes_vu .* nlp.active_impulse_mask;

% Compute the total delta-v in meters per second.
dv_total_mps = sum(DV_mag_active_vu) * nlp.VU * 1000;

% Compute the total time of flight in days.
tof_days = dec.t_nodes(end) * setup.constants.TU_days;

% Evaluate the constraint residuals for reporting.
[c_val, ceq_val] = constr_vleo_rpo_regularized_2d(y, nlp);

% Compute the minimum Moon altitude along the propagated trajectory.
r_moon_all_du = hypot(X_traj(:, 1) - 1 + nlp.mu, X_traj(:, 2));
moon_altitude_all_km = r_moon_all_du * nlp.DU - nlp.R_m;

% Read the circular departure radius and speed in physical units.
departure_radius_km = nlp.R_e + dec.departure_altitude_km;
departure_circular_velocity_kms = sqrt(nlp.muE / departure_radius_km);

% Read the departure pre- and post-impulse inertial speeds.
departure_pre_speed_kms = norm(x_departure_pre_in(3:4)) * nlp.VU;
departure_post_speed_kms = norm(x_departure_post_in(3:4)) * nlp.VU;

% Read the departure phase in degrees for the terminal logs.
departure_phase_deg = rad2deg(dec.phase_departure);

% Read the arrival Moon-relative distance and altitude.
arrival_relative_position_km = nlp.DU * (x_arrival_target_in(1:2) - x_moon_final_in(1:2));
arrival_moon_distance_km = norm(arrival_relative_position_km);
arrival_moon_altitude_km = arrival_moon_distance_km - nlp.R_m;

% Read the arrival Moon-relative pre- and post-insertion speeds.
arrival_pre_relative_velocity_kms = nlp.VU * (x_arrival_pre_in(3:4) - x_moon_final_in(3:4));
arrival_target_relative_velocity_kms = nlp.VU * (x_arrival_target_in(3:4) - x_moon_final_in(3:4));
arrival_pre_relative_speed_kms = norm(arrival_pre_relative_velocity_kms);
arrival_target_relative_speed_kms = norm(arrival_target_relative_velocity_kms);

% Store the compact departure diagnostics block.
departure_summary = struct;
departure_summary.leo_altitude_km = dec.departure_altitude_km;
departure_summary.leo_radius_km = departure_radius_km;
departure_summary.circular_velocity_kms = departure_circular_velocity_kms;
departure_summary.pre_impulse_speed_kms = departure_pre_speed_kms;
departure_summary.post_impulse_speed_kms = departure_post_speed_kms;
departure_summary.impulse_magnitude_mps = norm(dv_departure) * nlp.VU * 1000;
departure_summary.phase_deg = departure_phase_deg;

% Store the compact arrival diagnostics block.
arrival_summary = struct;
arrival_summary.final_moon_distance_km = arrival_moon_distance_km;
arrival_summary.final_moon_altitude_km = arrival_moon_altitude_km;
arrival_summary.pre_insertion_speed_rel_moon_kms = arrival_pre_relative_speed_kms;
arrival_summary.target_rpo_speed_rel_moon_kms = arrival_target_relative_speed_kms;
arrival_summary.impulse_magnitude_mps = norm(dv_arrival) * nlp.VU * 1000;

% Package the decoded solution and summary metrics.
solution = struct;
solution.tag = tag;
solution.y = y;
solution.dec = dec;
solution.X_nodes = X_nodes;
solution.x_departure_pre = x_departure_pre;
solution.x_rpo_target = x_rpo_target;
solution.t_traj = t_traj;
solution.X_traj = X_traj;
solution.DV_nodes_rot = DV_nodes_rot;
solution.DV_mag_nodes_vu = DV_mag_nodes_vu;
solution.DV_mag_active_vu = DV_mag_active_vu;
solution.dv_total_mps = dv_total_mps;
solution.tof_days = tof_days;
solution.max_abs_ceq = max(abs(ceq_val));
solution.max_c = max(c_val);
solution.min_moon_altitude_km = min(moon_altitude_all_km);
solution.theta_sun0 = dec.theta_sun0;
solution.phase_departure = dec.phase_departure;
solution.phase_rpo = dec.phase_rpo;
solution.departure_altitude_km = dec.departure_altitude_km;
solution.departure_summary = departure_summary;
solution.arrival_summary = arrival_summary;

end
