function [x_li_nd, gamma_li_nd] = libration_point_collinear_cr3bp(mu_nd, li_id)
%LIBRATION_POINT_COLLINEAR_CR3BP Locate a collinear libration point and its primary-relative distance.
%
% PROTOTYPE:
%   [x_li_nd, gamma_li_nd] = libration_point_collinear_cr3bp(mu_nd, li_id)
%
% DESCRIPTION:
%   Bracket and solve the collinear CR3BP equilibrium equation
%   in the rotating barycentric frame, with primaries at x=-mu_nd
%   and x=1-mu_nd. The second output is the Richardson length scale.
%
% INPUT:
%   mu_nd  Scalar secondary/total mass ratio in (0,0.5). [-]
%   li_id  Text identifying L1, L2 or L3 (case-insensitive).
%
% OUTPUT:
%   x_li_nd      Scalar libration-point x coordinate [DU].
%   gamma_li_nd  Positive distance to the secondary for L1/L2, or
%                  to the primary for L3 [DU].

    validateattributes(mu_nd, {'numeric'}, {'real', 'finite', 'scalar', '>', 0, '<', 0.5}, ...
        mfilename, 'mu_nd', 1);

    li_id = upper(string(li_id));
    valid_li = any(li_id == ["L1", "L2", "L3"]);
    if ~valid_li
        error('%s:InvalidLibrationPoint', mfilename, ...
            'li_id must be one of ''L1'', ''L2'', ''L3''.');
    end

    % Build the equilibrium residual and stay away from primary singularities.
    f_col = @(x) x ...
        - (1 - mu_nd) * (x + mu_nd) ./ abs(x + mu_nd).^3 ...
        - mu_nd * (x - 1 + mu_nd) ./ abs(x - 1 + mu_nd).^3;

    eps_sing = 1e-8;

    % Choose the appropriate interval relative to the two primaries.
    switch li_id
        case "L1"
            a = -mu_nd + eps_sing;
            b = 1 - mu_nd - eps_sing;
        case "L2"
            a = 1 - mu_nd + eps_sing;
            b = a + 5.0;
        otherwise % L3
            a = -5.0;
            b = -mu_nd - eps_sing;
    end

    % Find a finite sign-changing bracket before solving for the selected point.
    [ok_bracket, a_bracket, b_bracket] = local_find_bracket(f_col, a, b, 2000);
    if ~ok_bracket
        error('%s:BracketFailure', mfilename, ...
            'Unable to bracket a root for %s in [%.6f, %.6f].', li_id, a, b);
    end

    x_li_nd = fzero(f_col, [a_bracket, b_bracket]);

    switch li_id
        case "L1"
            gamma_li_nd = (1 - mu_nd) - x_li_nd;
        case "L2"
            gamma_li_nd = x_li_nd - (1 - mu_nd);
        otherwise % L3
            gamma_li_nd = -mu_nd - x_li_nd;
    end

    if ~(isfinite(gamma_li_nd) && gamma_li_nd > 0)
        error('%s:InvalidGamma', mfilename, ...
            'Computed gamma is invalid for %s (gamma=%.6e).', li_id, gamma_li_nd);
    end
end

function [ok, a_sel, b_sel] = local_find_bracket(fun_handle, a, b, n_samples)
%LOCAL_FIND_BRACKET Scan n_samples points in [a,b] for the first finite root bracket.
% fun_handle is vectorized; return ok=false and NaN endpoints if none is found.

    ok = false;
    a_sel = NaN;
    b_sel = NaN;

    x = linspace(a, b, n_samples);
    f = fun_handle(x);

    for k = 1:(n_samples - 1)
        f1 = f(k);
        f2 = f(k + 1);

        if ~isfinite(f1) || ~isfinite(f2)
            continue;
        end

        if f1 == 0
            ok = true;
            a_sel = x(k);
            b_sel = x(min(k + 1, n_samples));
            return;
        end

        if sign(f1) ~= sign(f2)
            ok = true;
            a_sel = x(k);
            b_sel = x(k + 1);
            return;
        end
    end
end
