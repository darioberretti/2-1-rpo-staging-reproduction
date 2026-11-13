function seed = richardson_seed_halo_se3bp(mu_nd, li_id, Az_target_nd, varargin)
%RICHARDSON_SEED_HALO_SE3BP Evaluate a third-order Richardson halo seed about L1 or L2.
%
% PROTOTYPE:
%   seed = richardson_seed_halo_se3bp(mu_nd, li_id, Az_target_nd, Name, Value, ...)
%
% DESCRIPTION:
%   Determine the local planar amplitude from the prescribed
%   vertical amplitude, evaluate the third-order series and transform
%   the local state to rotating barycentric coordinates. This is an
%   analytical seed, not a differentially corrected periodic orbit.
%
% INPUT:
%   mu_nd         Scalar secondary/total mass ratio in (0,0.5). [-]
%   li_id         L1 or L2 text, case-insensitive.
%   Az_target_nd  Positive target vertical amplitude [DU_SE].
%   Name, Value   Optional pairs:
%       north_south  Nonzero scalar; its sign selects north (+) or
%                      south (-), default +1.
%       t_nd         Evaluation epoch [TU_SE], default 0.
%       phase_rad    Series phase [rad]; default pi for L1 or 0 for L2.
%
% OUTPUT:
%   seed  Struct with x0_nd [6 x 1] = [x;y;z;vx;vy;vz] in the SE
%           rotating barycentric frame, positions [DU_SE], velocity
%           [VU_SE]; T_est_nd is the estimated period [TU_SE].
%           state0_li_nd uses the libration-centered gamma length
%           scale; Ax_li_nd/Az_li_nd use that same local scale.
%           Ax_est_nd/Az_target_nd use DU_SE. Also returns coeff
%           (richardson_coefficients_cr3bp), geometry, branch sign,
%           phase/epoch and frequency correction nu [-].

    validateattributes(mu_nd, {'numeric'}, {'real', 'finite', 'scalar', '>', 0, '<', 0.5}, ...
        mfilename, 'mu_nd', 1);

    li_id = upper(string(li_id));
    if ~any(li_id == ["L1", "L2"])
        error('%s:UnsupportedLibrationPoint', mfilename, ...
            'Only L1 and L2 are supported for this halo seed function.');
    end

    validateattributes(Az_target_nd, {'numeric'}, {'real', 'finite', 'scalar', '>', 0}, ...
        mfilename, 'Az_target_nd', 3);

    parser = inputParser;
    parser.addParameter('north_south', 1, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && abs(v) > 0);
    parser.addParameter('t_nd', 0, @(v) isnumeric(v) && isscalar(v) && isfinite(v));
    parser.addParameter('phase_rad', [], @(v) isempty(v) || (isnumeric(v) && isscalar(v) && isfinite(v)));
    parser.parse(varargin{:});
    opts = parser.Results;

    dm = sign(opts.north_south);
    t_nd = opts.t_nd;

    if isempty(opts.phase_rad)
        if li_id == "L1"
            phase_rad = pi;
        else
            phase_rad = 0;
        end
    else
        phase_rad = opts.phase_rad;
    end

    % Recover local amplitudes and the corrected series frequency.
    coeff = richardson_coefficients_cr3bp(mu_nd, li_id);

    Az_li_nd = Az_target_nd / coeff.gamma_li_nd;
    radicand = -(coeff.l2 * Az_li_nd^2 + coeff.Delta) / coeff.l1;
    if ~(isfinite(radicand) && radicand > 0)
        error('%s:InvalidAmplitude', mfilename, ...
            ['Invalid Ax_li radicand (%.6e). ', ...
             'Chosen Az may be outside third-order halo validity for %s.'], radicand, li_id);
    end

    Ax_li_nd = sqrt(radicand);

    nu = 1 + coeff.s1 * Ax_li_nd^2 + coeff.s2 * Az_li_nd^2;
    T_est_nd = 2 * pi / (coeff.omega_p * nu);

    tau1 = coeff.omega_p * nu * t_nd + phase_rad;
    dtau1dt = coeff.omega_p * nu;

    % Evaluate third-order positions and their analytic time derivatives.
    state_li_nd = zeros(6, 1);
    state_li_nd(1) = coeff.a21 * Ax_li_nd^2 + coeff.a22 * Az_li_nd^2 ...
        - Ax_li_nd * cos(tau1) ...
        + (coeff.a23 * Ax_li_nd^2 - coeff.a24 * Az_li_nd^2) * cos(2 * tau1) ...
        + (coeff.a31 * Ax_li_nd^3 - coeff.a32 * Ax_li_nd * Az_li_nd^2) * cos(3 * tau1);

    state_li_nd(2) = coeff.kappa * Ax_li_nd * sin(tau1) ...
        + (coeff.b21 * Ax_li_nd^2 - coeff.b22 * Az_li_nd^2) * sin(2 * tau1) ...
        + (coeff.b31 * Ax_li_nd^3 - coeff.b32 * Ax_li_nd * Az_li_nd^2) * sin(3 * tau1);

    state_li_nd(3) = dm * Az_li_nd * cos(tau1) ...
        + dm * coeff.d21 * Ax_li_nd * Az_li_nd * (cos(2 * tau1) - 3) ...
        + dm * (coeff.d32 * Az_li_nd * Ax_li_nd^2 - coeff.d31 * Az_li_nd^3) * cos(3 * tau1);

    state_li_nd(4) = dtau1dt * ( ...
        Ax_li_nd * sin(tau1) ...
        - 2 * (coeff.a23 * Ax_li_nd^2 - coeff.a24 * Az_li_nd^2) * sin(2 * tau1) ...
        - 3 * (coeff.a31 * Ax_li_nd^3 - coeff.a32 * Ax_li_nd * Az_li_nd^2) * sin(3 * tau1));

    state_li_nd(5) = dtau1dt * ( ...
        coeff.kappa * Ax_li_nd * cos(tau1) ...
        + 2 * (coeff.b21 * Ax_li_nd^2 - coeff.b22 * Az_li_nd^2) * cos(2 * tau1) ...
        + 3 * (coeff.b31 * Ax_li_nd^3 - coeff.b32 * Ax_li_nd * Az_li_nd^2) * cos(3 * tau1));

    state_li_nd(6) = dtau1dt * ( ...
        -dm * Az_li_nd * sin(tau1) ...
        - 2 * dm * coeff.d21 * Ax_li_nd * Az_li_nd * sin(2 * tau1) ...
        - 3 * dm * (coeff.d32 * Az_li_nd * Ax_li_nd^2 - coeff.d31 * Az_li_nd^3) * sin(3 * tau1));

    state_li_nd(abs(state_li_nd) < eps) = 0;

    % Map the local seed to barycentric units and retain the series diagnostics.
    state0_nd = local_li_to_synodic(li_id, coeff.gamma_li_nd, mu_nd, state_li_nd);

    seed = struct();
    seed.mu_nd = mu_nd;
    seed.li_id = char(li_id);
    seed.gamma_li_nd = coeff.gamma_li_nd;
    seed.x_li_nd = coeff.x_li_nd;

    seed.north_south = dm;
    seed.phase_rad = phase_rad;
    seed.t_nd = t_nd;

    seed.Az_target_nd = Az_target_nd;
    seed.Az_li_nd = Az_li_nd;
    seed.Ax_li_nd = Ax_li_nd;
    seed.Ax_est_nd = Ax_li_nd * coeff.gamma_li_nd;
    seed.nu = nu;
    seed.T_est_nd = T_est_nd;

    seed.state0_li_nd = state_li_nd;
    seed.x0_nd = state0_nd;

    seed.coeff = coeff;
end

function state0_nd = local_li_to_synodic(li_id, gamma_li_nd, mu_nd, state_li_nd)
%LOCAL_LI_TO_SYNODIC Scale an L1/L2-centered [6 x 1] state by gamma and shift its x origin.
% Input local lengths use gamma*DU; output uses barycentric DU and DU/TU.

    state0_nd = zeros(6, 1);

    switch li_id
        case "L1"
            state0_nd(1) = state_li_nd(1) * gamma_li_nd - mu_nd + 1 - gamma_li_nd;
        case "L2"
            state0_nd(1) = state_li_nd(1) * gamma_li_nd - mu_nd + 1 + gamma_li_nd;
        otherwise
            error('richardson_seed_halo_se3bp:UnsupportedLibrationPoint', ...
                'Only L1/L2 mapping is supported in this function.');
    end

    state0_nd(2:6) = gamma_li_nd * state_li_nd(2:6);
end
