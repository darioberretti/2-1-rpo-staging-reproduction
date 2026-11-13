function [records, diagnostics] = compute_seed(row, ctx)
%COMPUTE_SEED Propagate one plus-branch SEL2 seed and return its lunar entries.
%
% PROTOTYPE:
%   [records, diagnostics] = sel2moon.compute_seed(row, context)
%
% DESCRIPTION:
%   Perturb velocity along the transported unstable direction, propagate in
%   SE-CR3BP, map to EMB, then propagate in BCR4BP. Record all decreasing
%   Moon-SOI crossings. A stage-1 guard advances to stage 2 from its endpoint.
%   Dense integration samples are temporary and are not returned.
%
% INPUT:
%   row          Scalar phase index in 1:settings.n_seed.
%   ctx          Output of sel2moon.prepare_seeds. Each model uses its own
%                rotating [x,y,z,vx,vy,vz] state, DU, VU and TU.
%
% OUTPUT:
%   records      [H x 1 struct] Encounter rows for this seed, using the
%                22-column SEL2MoonDatabase schema. Relative time includes
%                both propagation stages; absolute time also includes halo
%                phase waiting. EM normalized times and physical days are
%                both stored. vInf* fields are local relative velocities
%                at the SOI boundary [km/s], not asymptotic excess velocities.
%   diagnostics  Compact encounter count and terminal-event IDs for each
%                stage (zero means the horizon was reached), plus whether
%                stage 1 stopped before its requested horizon.

cfg = ctx.cfg;
k_seed = row;
branch_name = "plus";

% Transport unstable direction from the monodromy eigenvector to current phase.
x_tau = ctx.x_hist(ctx.idx_tau(k_seed), :).';
Phi_tau = ctx.Phi_hist(:, :, ctx.idx_tau(k_seed));
v_tau = local_enforce_real(Phi_tau * ctx.v_u0);
v_tau = local_normalize(v_tau);

% Build one perturbed seed in SE coordinates.
x_seed_se = x_tau;
x_seed_se(4:6) = x_seed_se(4:6) + cfg.epsilon_vel_nd * v_tau(4:6);
theta_seed = ctx.theta_sun_seed(k_seed);

% Propagate the stage-1 segment in SE3BP with guard events.
ev_se = @(t, x) local_event_se_stop( ...
    t, x, theta_seed, cfg.ws_sun_em_nd, ctx.t_ratio_se_to_em, ...
    ctx.mu_se, cfg.mu_em_nd, ctx.r_earth_guard_se_nd, ...
    ctx.r_moon_guard_se_nd, ctx.r_outer_guard_se_nd, ...
    cfg.alpha_plus_theta_const, ctx.L_se_km, cfg.du_em_km, ctx.omega_se, cfg.omega_em_rad_s);
opts_se_stage1 = odeset(ctx.opts_se, 'Events', ev_se);
[t_first, xx_first_se, te_first, ~, ie_first] = ode113( ...
    @(t, x) rhs_SE3BP_3D(t, x, ctx.mu_se), ctx.t_eval_first, x_seed_se, opts_se_stage1);

require_complete(t_first, xx_first_se, te_first, ie_first, ...
    ctx.t_eval_first(end), [1, 2, 3], row, 1);

% Convert final stage-1 state from SE frame to EMB frame.
x_switch_se = xx_first_se(end, :).';
t_first_em_end = (t_first(end) * ctx.TU_se_days) / cfg.tu_em_days;
theta_switch = theta_seed + cfg.ws_sun_em_nd * t_first_em_end;
x_seed_emb = se_cr3bp_to_emb_bcr4bp_state( ...
    x_switch_se.', theta_switch, ...
    'alpha_plus_theta_sun_const_rad', cfg.alpha_plus_theta_const, ...
    'l_se_km', ctx.L_se_km, ...
    'l_em_km', cfg.du_em_km, ...
    'omega_se_rad_s', ctx.omega_se, ...
    'omega_em_rad_s', cfg.omega_em_rad_s, ...
    'mu_se_nd', ctx.mu_se, 'mu_em_nd', cfg.mu_em_nd, ...
    'se_secondary_center', 'emb');
x_seed_emb = x_seed_emb(:);

% Propagate stage-2 segment in BCR4BP and detect Moon-SOI entries.
rhs_bcr = @(t, x) rhs_BCR4BP_3D(t, x, cfg.mu_em_nd, theta_switch);
ev_bcr = @(t, x) local_event_bcr4bp_stop( ...
    t, x, cfg.mu_em_nd, cfg.r_moon_soi_nd, cfg.r_earth_guard_nd, ...
    cfg.r_outer_guard_nd, cfg.r_moon_guard_nd);
opts_bcr = odeset('RelTol', cfg.rel_tol, 'AbsTol', cfg.abs_tol, 'Events', ev_bcr);
[t_prop, xx_prop, te, ye, ie] = ode113(rhs_bcr, ctx.t_eval_prop, x_seed_emb, opts_bcr);

require_complete(t_prop, xx_prop, te, ie, ctx.t_eval_prop(end), [2,3,4], row, 2);
idx_moon_evt = find(ie == 1);

% Build row-level encounter records for each Moon-SOI entry.
n_enc = numel(idx_moon_evt);
enc_local = repmat(local_empty_enc_record(), n_enc, 1);
for j = 1:n_enc
    k_evt = idx_moon_evt(j);
    t_enc_rel_nd = t_first_em_end + te(k_evt);
    x_enc_emb_nd = ye(k_evt, 1:6);
    t_enc_abs_nd = ctx.seed_elapsed_em_nd(k_seed) + t_enc_rel_nd;
    m = local_compute_encounter_metrics(x_enc_emb_nd, t_enc_abs_nd, ...
        cfg.mu_em_nd, cfg.du_em_km, cfg.tu_em_days);

    rec = local_empty_enc_record();
    rec.seed_global_idx = row;
    rec.seed_phase_idx = k_seed;
    rec.branch = branch_name;
    rec.tau_seed_nd = ctx.tau_seed_nd(k_seed);
    rec.seed_time_from_theta0_nd = ctx.seed_elapsed_em_nd(k_seed);
    rec.theta_sun_seed_rad = theta_seed;
    rec.encounter_idx_seed = j;
    rec.n_encounters_seed = n_enc;
    rec.t_enc_rel_nd = t_enc_rel_nd;
    rec.t_enc_rel_days_em = t_enc_rel_nd * cfg.tu_em_days;
    rec.t_enc_abs_nd = t_enc_abs_nd;
    rec.t_enc_abs_days_em = t_enc_abs_nd * cfg.tu_em_days;
    rec.x_seed_se_nd = x_seed_se.';
    rec.x_seed_emb_nd = x_seed_emb.';
    rec.x_enc_emb_nd = x_enc_emb_nd;
    rec.v_sc_eci_km_s = m.v_sc_eci_km_s;
    rec.v_moon_eci_km_s = m.v_moon_eci_km_s;
    rec.vInfMoon_ECI_km_s = m.v_inf_eci_km_s;
    rec.vInfMoon_RTH_km_s = m.v_inf_rth_km_s;
    rec.vInfMagMoon_km_s = m.v_inf_mag_km_s;
    rec.thetaMoon_deg = m.theta_moon_deg;
    enc_local(j, 1) = rec;
end

% Retain only the encounter records and compact termination diagnostics.
records = enc_local;
first_event = 0;
if t_first(end) < ctx.t_eval_first(end)-1e-11, first_event = ie_first(end); end
second_event = 0;
if t_prop(end) < ctx.t_eval_prop(end)-1e-11, second_event = ie(end); end
diagnostics = struct('n_encounters', n_enc, ...
    'stage1_ended_early', t_first(end) < ctx.t_eval_first(end)-1e-11, ...
    'stage1_terminal_event', first_event, 'stage2_terminal_event', second_event);
end

function require_complete(time, states, event_time, event_ids, horizon, terminal_ids, row, stage)
%REQUIRE_COMPLETE Reject solver truncation without a valid terminal guard.
% time/states are the solver samples [TU, DU/VU]; event_time/event_ids are
% its event arrays. horizon is the requested final time [TU]; terminal_ids
% lists valid guard IDs. row/stage identify a failure in the error message.
terminal = ~isempty(event_ids) && ismember(event_ids(end), terminal_ids) ...
    && abs(time(end)-event_time(end)) <= 1e-10*max(1,abs(horizon));
if any(~isfinite(states),'all') || (time(end) ~= horizon && ~terminal)
    error('sel2moon:IncompletePropagation', ...
        'Seed %d stage %d stopped at %.17g before horizon %.17g without a terminal guard.', ...
        row, stage, time(end), horizon);
end
end

function m = local_compute_encounter_metrics(x_enc_emb_nd, t_abs_nd, mu_em_nd, du_em_km, tu_em_days)
%LOCAL_COMPUTE_ENCOUNTER_METRICS Express boundary velocities in ECI and lunar RTH.
%
% PROTOTYPE:
%   m = local_compute_encounter_metrics(x_enc_emb_nd, t_abs_nd, ...
%       mu_em_nd, du_em_km, tu_em_days)
%
% DESCRIPTION:
%   Subtract the Moon's Earth-centered inertial (ECI) velocity from the
%   spacecraft velocity and project onto lunar radial/transverse/normal
%   axes (RTH). These are local SOI-boundary velocities.
%
% INPUT:
%   x_enc_emb_nd  [1 x 6] EMB rotating [x,y,z,vx,vy,vz] in EM DU/VU.
%   t_abs_nd      Absolute encounter time [TU_EM], including halo waiting.
%   mu_em_nd      Earth-Moon mass parameter [-].
%   du_em_km      Distance unit [km].
%   tu_em_days    Time unit [days].
%
% OUTPUT:
%   m             Struct of [1 x 3] spacecraft/Moon ECI velocities and
%                 relative ECI/RTH velocities [km/s], their scalar norm
%                 [km/s], and wrapped lunar phase theta_moon_deg [deg].
vu_em_km_s = du_em_km / (tu_em_days * 86400);

sc_inert_2d = rot2in([x_enc_emb_nd(1), x_enc_emb_nd(2), ...
    x_enc_emb_nd(4), x_enc_emb_nd(5)], t_abs_nd, mu_em_nd);
v_sc_eci_km_s = [sc_inert_2d(3), sc_inert_2d(4), x_enc_emb_nd(6)] * vu_em_km_s;

moon_rot_2d = [1 - mu_em_nd, 0, 0, 0];
moon_inert_2d = rot2in(moon_rot_2d, t_abs_nd, mu_em_nd);
v_moon_eci_km_s = [moon_inert_2d(3), moon_inert_2d(4), 0] * vu_em_km_s;
r_moon_eci_km = [moon_inert_2d(1), moon_inert_2d(2), 0] * du_em_km;

v_inf_eci_km_s = v_sc_eci_km_s - v_moon_eci_km_s;

r_hat = r_moon_eci_km / max(norm(r_moon_eci_km), realmin);
h_hat = [0, 0, 1];
t_hat = cross(h_hat, r_hat);
t_hat = t_hat / max(norm(t_hat), realmin);

v_inf_rth_km_s = [dot(v_inf_eci_km_s, r_hat), dot(v_inf_eci_km_s, t_hat), ...
    dot(v_inf_eci_km_s, h_hat)];

m = struct();
m.v_sc_eci_km_s = v_sc_eci_km_s;
m.v_moon_eci_km_s = v_moon_eci_km_s;
m.v_inf_eci_km_s = v_inf_eci_km_s;
m.v_inf_rth_km_s = v_inf_rth_km_s;
m.v_inf_mag_km_s = norm(v_inf_rth_km_s);
m.theta_moon_deg = rad2deg(local_wrap_to_2pi(t_abs_nd));
end

function v = local_enforce_real(v_in)
%LOCAL_ENFORCE_REAL Reject a significantly complex unstable direction.
if norm(imag(v_in))/max(norm(v_in),realmin) < 1e-8
    v = real(v_in);
else
    error('sel2moon:UnstableDirection','The transported unstable direction is not real.');
end
end

function v = local_normalize(v_in)
%LOCAL_NORMALIZE Normalize the full six-component state direction.
nv = norm(v_in);
if nv <= eps
    error('sel2moon:UnstableDirection','Cannot normalize a near-zero direction.');
end
v = v_in/nv;
end

function y = local_wrap_to_2pi(x)
%LOCAL_WRAP_TO_2PI Wrap an angle in radians to [0,2*pi).
y = mod(x, 2*pi);
y(y < 0) = y(y < 0) + 2*pi;
end

function [value, isterminal, direction] = local_event_se_stop( ...
    t, xx, theta_seed, ws_sun_em_nd, t_ratio_se_to_em, ...
    mu_se_nd, mu_em_nd, r_earth_guard_se_nd, r_moon_guard_se_nd, r_outer_guard_se_nd, ...
    alpha_plus_theta_const, l_se_km, l_em_km, omega_se_rad_s, omega_em_rad_s)
%LOCAL_EVENT_SE_STOP Detect stage-1 Earth, Moon and outer-radius guards.
% The [6 x 1] state xx uses Sun-Earth rotating DU/VU; t uses TU_SE.
% theta_seed and alpha_plus_theta_const are angles [rad]; ws_sun_em_nd
% is the Sun's rate [rad/TU_EM]; t_ratio_se_to_em converts TU_SE to TU_EM.
% Mass parameters are dimensionless, radii use DU_SE, lengths use km,
% and angular rates use rad/s. The explicit EMB mapping locates the Moon.
% Return three signed distances [DU_SE], terminal flags and negative
% crossing directions in Earth/Moon/outer order.
earth_center_se = [1 - mu_se_nd, 0, 0];
theta_abs = theta_seed + ws_sun_em_nd * (t * t_ratio_se_to_em);
moon_state_se = emb_bcr4bp_to_se_cr3bp_state( ...
    [1 - mu_em_nd, 0, 0, 0, 0, 0], theta_abs, ...
    'alpha_plus_theta_sun_const_rad', alpha_plus_theta_const, ...
    'l_se_km', l_se_km, ...
    'l_em_km', l_em_km, ...
    'omega_se_rad_s', omega_se_rad_s, ...
    'omega_em_rad_s', omega_em_rad_s, ...
    'mu_se_nd', mu_se_nd, 'mu_em_nd', mu_em_nd, ...
    'se_secondary_center', 'emb');
moon_center_se = moon_state_se(1, 1:3);

r_earth = norm(xx(1:3).' - earth_center_se);
r_moon = norm(xx(1:3).' - moon_center_se);
r_norm = norm(xx(1:3));

value = [ ...
    r_earth - r_earth_guard_se_nd; ...
    r_moon - r_moon_guard_se_nd; ...
    r_outer_guard_se_nd - r_norm];
isterminal = [1; 1; 1];
direction = [-1; -1; -1];
end

function [value, isterminal, direction] = local_event_bcr4bp_stop( ...
    ~, xx, mu_em_nd, r_moon_soi_nd, r_earth_guard_nd, r_outer_guard_nd, r_moon_guard_nd)
%LOCAL_EVENT_BCR4BP_STOP Log lunar entries and detect stage-2 guards.
% xx is a [6 x 1] EMB rotating state in EM DU/VU; all radii use DU_EM.
% Return four signed distances [DU_EM], terminal flags and negative
% crossing directions in SOI/Earth/outer/Moon order. SOI is nonterminal.
r_earth = norm(xx(1:3).' - [-mu_em_nd, 0, 0]);
r_moon = norm(xx(1:3).' - [1 - mu_em_nd, 0, 0]);
r_norm = norm(xx(1:3));

value = [ ...
    r_moon - r_moon_soi_nd; ...
    r_earth - r_earth_guard_nd; ...
    r_outer_guard_nd - r_norm; ...
    r_moon - r_moon_guard_nd];
isterminal = [0; 1; 1; 1];
direction = [-1; -1; -1; -1];
end

function rec = local_empty_enc_record()
%LOCAL_EMPTY_ENC_RECORD Define one scalar row of the 22-column encounter table.
rec = struct( ...
    'encounter_global_idx', 0, ...
    'seed_global_idx', 0, ...
    'seed_phase_idx', 0, ...
    'branch', "", ...
    'tau_seed_nd', 0, ...
    'seed_time_from_theta0_nd', 0, ...
    'theta_sun_seed_rad', 0, ...
    'encounter_idx_seed', 0, ...
    'n_encounters_seed', 0, ...
    't_enc_rel_nd', 0, ...
    't_enc_rel_days_em', 0, ...
    't_enc_abs_nd', 0, ...
    't_enc_abs_days_em', 0, ...
    'x_seed_se_nd', zeros(1, 6), ...
    'x_seed_emb_nd', zeros(1, 6), ...
    'x_enc_emb_nd', zeros(1, 6), ...
    'v_sc_eci_km_s', zeros(1, 3), ...
    'v_moon_eci_km_s', zeros(1, 3), ...
    'vInfMoon_ECI_km_s', zeros(1, 3), ...
    'vInfMoon_RTH_km_s', zeros(1, 3), ...
    'vInfMagMoon_km_s', 0, ...
    'thetaMoon_deg', 0);
end
