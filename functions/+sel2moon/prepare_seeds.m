function context = prepare_seeds(settings)
%PREPARE_SEEDS Reconstruct the halo and transport its unstable direction.
%
% PROTOTYPE:
%   context = sel2moon.prepare_seeds(settings)
%
% DESCRIPTION:
%   Correct the SEL2 halo, integrate its state-transition matrix over one
%   period, and select the largest-magnitude unstable multiplier. Align the
%   real eigenvector with the supplied reference to fix the plus-branch sign.
%   Seed phases exclude the period endpoint; Sun phase includes halo elapsed
%   time. The STM maps initial Sun-Earth rotating-state perturbations to
%   perturbations at each phase in the same normalization.
%
% INPUT:
%   settings  Frozen profile from build_sel2_to_moon_database.m. halo contains
%             the correction/model settings; reference_direction is [6 x 1].
%
% OUTPUT:
%   context   Orbit, monodromy, sampled SE states [N+1 x 6], STM [6 x 6 x N+1],
%             phase/epoch arrays [N x 1], and integration grids. States use
%             [x,y,z,vx,vy,vz] in the corresponding rotating DU/VU; N=n_seed.
%             SE (Sun-Earth) and EM (Earth-Moon) time units are distinct,
%             each equal to 1/omega.

% Reconstruct the orbit using the same configuration as the halo study.
cfg = struct('meta', struct('name', 'sel2_to_moon'), ...
    'sel2_to_rpo', struct('halo', settings.halo));
halo_result = sel2rpo.compute_halo(cfg);
halo = halo_result.data;
c = halo.constants;
period = halo.solution.period_nd;
initial = halo.solution.x0_nd(:);
time_unit_se_days = (1 / c.omega_SE_rad_s) / 86400;

% Include the full-period endpoint only for the monodromy calculation.
theta0 = mod(settings.theta_halo0_nd, period);
dtau = period / settings.n_seed;
seed_elapsed_se = (0:settings.n_seed-1).' * dtau;
seed_elapsed_em = seed_elapsed_se * time_unit_se_days / settings.tu_em_days;
phases = mod(theta0 + seed_elapsed_se, period);
sun_phases = settings.theta_sun0_rad + settings.ws_sun_em_nd * seed_elapsed_em;
time_grid = unique([0; phases; period], 'sorted');
options = odeset('RelTol', settings.rel_tol, 'AbsTol', settings.abs_tol);
identity = eye(6);
[time, extended] = ode113(@(t, x) rhs_SE3BP_3D_STM(t, x, c.mu_SE_nd), ...
    time_grid, [initial; identity(:)], options);
if time(end) ~= period || any(~isfinite(extended), 'all')
    error('sel2moon:IncompleteSTM', 'Halo STM propagation did not reach the full period.');
end
states = extended(:, 1:6);
transitions = reshape(extended(:, 7:42).', 6, 6, []);
monodromy = transitions(:, :, end);
[vectors, values] = eig(monodromy);
multipliers = diag(values);
[~, index] = max(abs(multipliers));
direction = vectors(:, index);
if norm(imag(direction)) / max(norm(direction), realmin) >= 1e-8
    error('sel2moon:UnstableDirection', 'The unstable direction is not real.');
end
direction = real(direction);
direction = direction / norm(direction);
if dot(direction, settings.reference_direction) < 0, direction = -direction; end
[present, phase_indices] = ismember(phases, time);
if ~all(present)
    error('sel2moon:PhaseGrid', 'Halo integration did not return all requested seed phases.');
end

% Keep normalized frame parameters and grids explicit for per-seed propagation.
context = struct('cfg', settings, 'mu_se', c.mu_SE_nd, 'L_se_km', c.L_SE_km, ...
    'omega_se', c.omega_SE_rad_s, 'TU_se_days', time_unit_se_days, ...
    'x_hist', states, 'Phi_hist', transitions, 'idx_tau', phase_indices, ...
    'v_u0', direction, 'tau_seed_nd', phases, 'seed_elapsed_em_nd', seed_elapsed_em, ...
    'theta_sun_seed', sun_phases, 'opts_se', options, ...
    't_ratio_se_to_em', time_unit_se_days / settings.tu_em_days, ...
    'r_earth_guard_se_nd', settings.r_earth_guard_nd * (settings.du_em_km / c.L_SE_km), ...
    'r_moon_guard_se_nd', settings.r_moon_guard_nd * (settings.du_em_km / c.L_SE_km), ...
    'r_outer_guard_se_nd', settings.r_outer_guard_nd, ...
    't_eval_first', linspace(0, settings.t_first_days / time_unit_se_days, ...
        settings.n_samples_first).', ...
    't_eval_prop', linspace(0, settings.t_second_days / settings.tu_em_days, ...
        settings.n_samples_second).', ...
    't_prop_second_nd', settings.t_second_days / settings.tu_em_days);
context.orbit = struct('x0_se_nd', initial, 'period_se_nd', period, ...
    'period_se_days', period * time_unit_se_days, 'mu_se_nd', c.mu_SE_nd, ...
    'l_se_km', c.L_SE_km, 'omega_se_rad_s', c.omega_SE_rad_s);
context.unstable = struct('monodromy_M', monodromy, 'lambda_u', multipliers(index), ...
    'eigvals_all', multipliers, 'unstable_index', index, 'direction', direction);
end
