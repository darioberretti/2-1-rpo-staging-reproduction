function coeff = richardson_coefficients_cr3bp(mu_nd, li_id)
%RICHARDSON_COEFFICIENTS_CR3BP Evaluate third-order Richardson coefficients for a collinear point.
%
% PROTOTYPE:
%   coeff = richardson_coefficients_cr3bp(mu_nd, li_id)
%
% DESCRIPTION:
%   Use the point location and gamma scale to evaluate the
%   third-order amplitude/frequency expansion. Formula reference:
%   Richardson (1980), Celestial Mechanics 22, 241-253,
%   DOI: 10.1007/BF01229511.
%
% INPUT:
%   mu_nd  Scalar secondary/total mass ratio in (0,0.5). [-]
%   li_id  L1, L2 or L3 text, case-insensitive.
%
% OUTPUT:
%   coeff  Struct with barycentric x_li_nd/gamma_li_nd [DU], mass
%            ratio and point ID; c1-c4, lambda_r/omega_p/omega_v, kappa,
%            Delta, d1/d2 and third-order a/b/d coefficients. s1/s2
%            correct frequency and l1/l2 couple amplitudes. Expansion
%            coefficients use local gamma-scaled coordinates and
%            normalized CR3BP time, as consumed by
%            richardson_seed_halo_se3bp.

    validateattributes(mu_nd, {'numeric'}, {'real', 'finite', 'scalar', '>', 0, '<', 0.5}, ...
        mfilename, 'mu_nd', 1);

    li_id = upper(string(li_id));
    if ~any(li_id == ["L1", "L2", "L3"])
        error('%s:InvalidLibrationPoint', mfilename, ...
            'li_id must be one of ''L1'', ''L2'', ''L3''.');
    end

    % Compute equilibrium geometry and local potential coefficients.
    [x_li_nd, gamma_li_nd] = libration_point_collinear_cr3bp(mu_nd, li_id);

    c1 = local_cn(mu_nd, li_id, gamma_li_nd, 1);
    c2 = local_cn(mu_nd, li_id, gamma_li_nd, 2);
    c3 = local_cn(mu_nd, li_id, gamma_li_nd, 3);
    c4 = local_cn(mu_nd, li_id, gamma_li_nd, 4);

    % Build linear frequencies and check denominators before nonlinear terms.
    lambda_r = sqrt(0.5 * (2 - c2 + sqrt(9 * c2^2 - 8 * c2)));
    omega_p = lambda_r;
    omega_v = sqrt(c2);

    kappa = (lambda_r^2 + 1 + 2 * c2) / (2 * lambda_r);
    Delta = lambda_r^2 - c2;

    d1 = 3 * lambda_r^2 / kappa * (kappa * (6 * lambda_r^2 - 1) - 2 * lambda_r);
    d2 = 8 * lambda_r^2 / kappa * (kappa * (11 * lambda_r^2 - 1) - 2 * lambda_r);

    local_require_nonzero(1 + 2 * c2, '1+2*c2');
    local_require_nonzero(d1, 'd1');
    local_require_nonzero(d2, 'd2');
    local_require_nonzero(lambda_r, 'lambda_r');

    % Evaluate second- and third-order position-series coefficients.
    a21 = (3 * c3 * (kappa^2 - 2)) / (4 * (1 + 2 * c2));
    a22 = (3 * c3) / (4 * (1 + 2 * c2));
    a23 = -(3 * c3 * lambda_r) / (4 * kappa * d1) * ...
        (3 * kappa^3 * lambda_r - 6 * kappa * (kappa - lambda_r) + 4);
    a24 = -(3 * c3 * lambda_r) / (4 * kappa * d1) * (2 + 3 * kappa * lambda_r);

    b21 = -(3 * c3 * lambda_r) / (2 * d1) * (3 * kappa * lambda_r - 4);
    b22 = (3 * c3 * lambda_r) / d1;

    d21 = -c3 / (2 * lambda_r^2);

    a31 = -(9 * lambda_r) / (4 * d2) * ...
        (4 * c3 * (kappa * a23 - b21) + kappa * c4 * (4 + kappa^2)) + ...
        (9 * lambda_r^2 + 1 - c2) / (2 * d2) * ...
        (3 * c3 * (2 * a23 - kappa * b21) + c4 * (2 + 3 * kappa^2));

    a32 = -(9 * lambda_r) / (4 * d2) * ...
        (4 * c3 * (kappa * a24 - b22) + kappa * c4) - ...
        3 / (2 * d2) * (9 * lambda_r^2 + 1 - c2) * ...
        (c3 * (kappa * b22 + d21 - 2 * a24) - c4);

    b31 = 3 / (8 * d2) * ( ...
        8 * lambda_r * (3 * c3 * (kappa * b21 - 2 * a23) - c4 * (2 + 3 * kappa^2)) + ...
        (9 * lambda_r^2 + 1 + 2 * c2) * (4 * c3 * (kappa * a23 - b21) + kappa * c4 * (4 + kappa^2)) ...
        );

    b32 = 1 / d2 * ( ...
        9 * lambda_r * (c3 * (kappa * b22 + d21 - 2 * a24) - c4) + ...
        3 / 8 * (9 * lambda_r^2 + 1 + 2 * c2) * (4 * c3 * (kappa * a24 - b22) + kappa * c4) ...
        );

    d31 = 3 / (64 * lambda_r^2) * (4 * c3 * a24 + c4);
    d32 = 3 / (64 * lambda_r^2) * (4 * c3 * (a23 - d21) + c4 * (4 + kappa^2));

    % Evaluate frequency corrections and the coupled amplitude relation.
    s_denom = 2 * lambda_r * (lambda_r * (1 + kappa^2) - 2 * kappa);
    local_require_nonzero(s_denom, 's_denom');

    s1 = (1 / s_denom) * ( ...
        3 / 2 * c3 * (2 * a21 * (kappa^2 - 2) - a23 * (kappa^2 + 2) - 2 * kappa * b21) - ...
        3 / 8 * c4 * (3 * kappa^4 - 8 * kappa^2 + 8) ...
        );

    s2 = (1 / s_denom) * ( ...
        3 / 2 * c3 * (2 * a22 * (kappa^2 - 2) + a24 * (kappa^2 + 2) + 2 * kappa * b22 + 5 * d21) + ...
        3 / 8 * c4 * (12 - kappa^2) ...
        );

    l1 = -3 / 2 * c3 * (2 * a21 + a23 + 5 * d21) - ...
        3 / 8 * c4 * (12 - kappa^2) + 2 * lambda_r^2 * s1;

    l2 = 3 / 2 * c3 * (a24 - 2 * a22) + 9 / 8 * c4 + 2 * lambda_r^2 * s2;

    % Return named coefficients alongside the geometry used to derive them.
    coeff = struct();
    coeff.mu_nd = mu_nd;
    coeff.li_id = char(li_id);
    coeff.x_li_nd = x_li_nd;
    coeff.gamma_li_nd = gamma_li_nd;

    coeff.c1 = c1;
    coeff.c2 = c2;
    coeff.c3 = c3;
    coeff.c4 = c4;

    coeff.lambda_r = lambda_r;
    coeff.omega_p = omega_p;
    coeff.omega_v = omega_v;
    coeff.kappa = kappa;
    coeff.Delta = Delta;
    coeff.d1 = d1;
    coeff.d2 = d2;

    coeff.a21 = a21;
    coeff.a22 = a22;
    coeff.a23 = a23;
    coeff.a24 = a24;
    coeff.b21 = b21;
    coeff.b22 = b22;
    coeff.d21 = d21;

    coeff.a31 = a31;
    coeff.a32 = a32;
    coeff.b31 = b31;
    coeff.b32 = b32;
    coeff.d31 = d31;
    coeff.d32 = d32;

    coeff.s1 = s1;
    coeff.s2 = s2;
    coeff.l1 = l1;
    coeff.l2 = l2;
end

function cn = local_cn(mu_nd, li_id, gamma_li_nd, n)
%LOCAL_CN Evaluate order-n local potential coefficient for the selected point.
% mu_nd is the mass ratio and gamma_li_nd is the primary-relative distance [DU].

    if li_id == "L1" && n == 1
        cn = (mu_nd - 1 + gamma_li_nd) / gamma_li_nd;
        return;
    end

    if li_id == "L2" && n == 1
        cn = (mu_nd - 1 - gamma_li_nd) / gamma_li_nd;
        return;
    end

    if li_id == "L3" && n == 1
        cn = (mu_nd + gamma_li_nd) / gamma_li_nd;
        return;
    end

    switch li_id
        case "L1"
            denom = (1 - gamma_li_nd)^(n + 1);
            cn = gamma_li_nd^(-3) * ( ...
                mu_nd + (-1)^n * (1 - mu_nd) * gamma_li_nd^(n + 1) / denom ...
                );
        case "L2"
            denom = (1 + gamma_li_nd)^(n + 1);
            cn = gamma_li_nd^(-3) * ( ...
                (-1)^n * mu_nd + (-1)^n * (1 - mu_nd) * gamma_li_nd^(n + 1) / denom ...
                );
        otherwise % L3
            denom = (1 + gamma_li_nd)^(n + 1);
            cn = gamma_li_nd^(-3) * ( ...
                1 - mu_nd + mu_nd * gamma_li_nd^(n + 1) / denom ...
                );
    end
end

function local_require_nonzero(value, label)
%LOCAL_REQUIRE_NONZERO Reject a nonfinite or near-zero scalar formula denominator (threshold 1e-14).

    if ~isfinite(value) || abs(value) < 1e-14
        error('richardson_coefficients_cr3bp:SingularCoefficient', ...
            'Coefficient denominator %s is singular (%.6e).', label, value);
    end
end
