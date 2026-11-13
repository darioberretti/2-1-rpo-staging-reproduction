function sol = diffcorr_halo_halfperiod_se3bp(x0_seed_nd, mu_se_nd, varargin)
%DIFFCORR_HALO_HALFPERIOD_SE3BP Correct a symmetric halo at fixed z0 using state/STM shooting.
%
% PROTOTYPE:
%   sol = diffcorr_halo_halfperiod_se3bp(x0_seed_nd, mu_se_nd, ...
%       'half_period_guess_nd', T_half, Name, Value, ...)
%
% DESCRIPTION:
%   Fix y0=vx0=vz0=0 and the seed z0. Solve for x0, vy0 and
%   half-period so [y,vx,vz] vanishes at the half-period. Build the
%   Newton Jacobian from an identity-initialized STM and terminal flow.
%   Report solver convergence separately from full-period closure.
%
% INPUT:
%   x0_seed_nd  [6 elements] Seed [x,y,z,vx,vy,vz] in the rotating
%                 SE barycentric frame, positions [DU_SE], velocity [VU_SE].
%   mu_se_nd    Scalar secondary/total mass ratio [-].
%   Name, Value Optional pairs (half_period_guess_nd is required):
%       half_period_guess_nd  Initial half-period [TU_SE], default [].
%       max_iter              Maximum Newton iterations, default 25.
%       tol_res_inf           Mixed normalized residual norm, default 1e-11.
%       tol_step_inf          Mixed normalized step norm, default 1e-12.
%       damping               Newton step factor in (0,1], default 1.
%       use_pinv              Use pinv below rcond threshold, default true.
%       cond_rcond_min        Jacobian rcond threshold, default 1e-12.
%       ode_rel_tol, ode_abs_tol  Integration tolerances, default 1e-12.
%       min_half_period_nd    Positive lower bound [TU_SE], default 1e-6.
%       verbose               Print iteration residuals, default false.
%
% OUTPUT:
%   sol  Struct with success/message; corrected x0_corrected_nd
%          [6 x 1], T_half_nd/T_nd [TU_SE], xf_half_nd [6 x 1], and
%          Phi_half [6 x 6] mapping initial to half-period state changes.
%          residual_final_nd [3 x 1] is [y;vx;vz]; closure_vec_nd [6 x 1]
%          is full-period endpoint minus initial state. Both have inf
%          norms in normalized components. trajectory contains t_nd
%          [M x 1, TU_SE] and xx_nd [M x 6]. Also retains seed, options,
%          iteration count, residual_history/step_history [3 x K] and
%          rcond_history [K x 1]. M is the adaptive ODE output count;
%          K is the number of performed Newton iterations.

    validateattributes(x0_seed_nd, {'numeric'}, {'real', 'finite', 'vector', 'numel', 6}, ...
        mfilename, 'x0_seed_nd', 1);
    validateattributes(mu_se_nd, {'numeric'}, {'real', 'finite', 'scalar', '>', 0, '<', 0.5}, ...
        mfilename, 'mu_se_nd', 2);

    parser = inputParser;
    parser.addParameter('half_period_guess_nd', [], ...
        @(v) isempty(v) || (isnumeric(v) && isscalar(v) && isfinite(v) && v > 0));
    parser.addParameter('max_iter', 25, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v >= 1 && mod(v,1) == 0);
    parser.addParameter('tol_res_inf', 1e-11, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v > 0);
    parser.addParameter('tol_step_inf', 1e-12, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v > 0);
    parser.addParameter('damping', 1.0, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v > 0 && v <= 1);
    parser.addParameter('use_pinv', true, @(v) islogical(v) || (isnumeric(v) && isscalar(v)));
    parser.addParameter('cond_rcond_min', 1e-12, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v > 0);
    parser.addParameter('ode_rel_tol', 1e-12, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v > 0);
    parser.addParameter('ode_abs_tol', 1e-12, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v > 0);
    parser.addParameter('min_half_period_nd', 1e-6, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v > 0);
    parser.addParameter('verbose', false, @(v) islogical(v) || (isnumeric(v) && isscalar(v)));
    parser.parse(varargin{:});
    opts = parser.Results;
    if isempty(opts.half_period_guess_nd)
        error('%s:MissingHalfPeriodGuess', mfilename, ...
            'Option ''half_period_guess_nd'' is required and must be positive.');
    end

    x0_seed_nd = x0_seed_nd(:);
    z0_fixed_nd = x0_seed_nd(3);

    % Initialize the free variables [x0,vy0,T_half] and iteration histories.
    u = [x0_seed_nd(1); x0_seed_nd(5); opts.half_period_guess_nd];
    u(3) = max(u(3), opts.min_half_period_nd);

    residual_history = zeros(3, opts.max_iter);
    step_history = zeros(3, opts.max_iter);
    rcond_history = zeros(opts.max_iter, 1);

    success = false;
    message = 'max_iter_reached';
    k_done = 0;

    % Propagate the trial state and correct its half-period symmetry residual.
    for k = 1:opts.max_iter
        k_done = k;
        x0_curr_nd = [u(1); 0; z0_fixed_nd; 0; u(2); 0];
        tf_nd = u(3);

        [ok_prop, prop_data, prop_msg] = local_propagate_halfperiod(x0_curr_nd, tf_nd, mu_se_nd, opts);
        if ~ok_prop
            success = false;
            message = prop_msg;
            break;
        end

        xf_half_nd = prop_data.xf_half_nd;
        Phi_half = prop_data.Phi_half;

        F = [xf_half_nd(2); xf_half_nd(4); xf_half_nd(6)];
        residual_history(:, k) = F;

        if opts.verbose
            fprintf('Iter %02d: ||F||_inf = %.3e\n', k, norm(F, inf));
        end

        if norm(F, inf) <= opts.tol_res_inf
            success = true;
            message = 'residual_tolerance_met';
            break;
        end

        % Differentiate symmetry conditions with respect to x0, vy0 and final time.
        f_tf = rhs_SE3BP_3D(tf_nd, xf_half_nd, mu_se_nd);
        J = [Phi_half([2, 4, 6], 1), Phi_half([2, 4, 6], 5), f_tf([2, 4, 6])];

        rc = rcond(J);
        rcond_history(k) = rc;

        if ~all(isfinite(J(:)))
            success = false;
            message = 'non_finite_jacobian';
            break;
        end

        if opts.use_pinv && rc < opts.cond_rcond_min
            dU = -pinv(J) * F;
        else
            dU = -J \ F;
        end

        if any(~isfinite(dU))
            success = false;
            message = 'invalid_newton_step';
            break;
        end

        dU = opts.damping * dU;
        step_history(:, k) = dU;

        u = u + dU;
        u(3) = max(u(3), opts.min_half_period_nd);

        if norm(dU, inf) <= opts.tol_step_inf
            if norm(F, inf) <= 10 * opts.tol_res_inf
                success = true;
                message = 'step_tolerance_met';
            else
                success = false;
                message = 'step_small_residual_not_met';
            end
            break;
        end
    end

    % Retain performed iterations and re-evaluate the final corrected half-period.
    residual_history = residual_history(:, 1:k_done);
    step_history = step_history(:, 1:k_done);
    rcond_history = rcond_history(1:k_done);

    x0_corrected_nd = [u(1); 0; z0_fixed_nd; 0; u(2); 0];
    T_half_nd = u(3);
    T_nd = 2 * T_half_nd;

    [ok_half, half_data, half_msg] = local_propagate_halfperiod(x0_corrected_nd, T_half_nd, mu_se_nd, opts);
    if ~ok_half
        success = false;
        message = half_msg;
        xf_half_nd = NaN(6, 1);
        Phi_half = NaN(6, 6);
        F_final = NaN(3, 1);
    else
        xf_half_nd = half_data.xf_half_nd;
        Phi_half = half_data.Phi_half;
        F_final = [xf_half_nd(2); xf_half_nd(4); xf_half_nd(6)];
    end

    rhs = @(t, x) rhs_SE3BP_3D(t, x, mu_se_nd);
    ode_opts = odeset('RelTol', opts.ode_rel_tol, 'AbsTol', opts.ode_abs_tol);

    % Independently integrate a full period to report the actual closure error.
    closure_vec_nd = NaN(6, 1);
    [tt_full_nd, xx_full_nd] = deal([]);
    try
        [tt_full_nd, xx_full_nd] = ode113(rhs, [0, T_nd], x0_corrected_nd, ode_opts);
        if ~isempty(xx_full_nd)
            closure_vec_nd = xx_full_nd(end, :).'-x0_corrected_nd;
        else
            success = false;
            if strcmp(message, 'residual_tolerance_met') || strcmp(message, 'step_tolerance_met')
                message = 'full_period_propagation_empty';
            end
        end
    catch
        success = false;
        if strcmp(message, 'residual_tolerance_met') || strcmp(message, 'step_tolerance_met')
            message = 'full_period_propagation_failed';
        end
    end

    % Keep solver status, geometric closure and full histories separately.
    sol = struct();
    sol.success = success;
    sol.message = message;

    sol.x0_seed_nd = x0_seed_nd;
    sol.x0_corrected_nd = x0_corrected_nd;
    sol.z0_fixed_nd = z0_fixed_nd;

    sol.T_half_nd = T_half_nd;
    sol.T_nd = T_nd;

    sol.residual_final_nd = F_final;
    sol.residual_inf_nd = norm(F_final, inf);
    sol.closure_vec_nd = closure_vec_nd;
    sol.closure_inf_nd = norm(closure_vec_nd, inf);

    sol.xf_half_nd = xf_half_nd;
    sol.Phi_half = Phi_half;

    sol.iterations = k_done;
    sol.residual_history = residual_history;
    sol.step_history = step_history;
    sol.rcond_history = rcond_history;

    sol.trajectory = struct('t_nd', tt_full_nd, 'xx_nd', xx_full_nd);
    sol.options = opts;
end

function [ok, data, msg] = local_propagate_halfperiod(x0_nd, tf_nd, mu_se_nd, opts)
%LOCAL_PROPAGATE_HALFPERIOD Integrate the SE state and identity-initialized STM to a trial half-period.
%
% PROTOTYPE:
%   [ok, data, msg] = local_propagate_halfperiod(x0_nd, tf_nd, mu_se_nd, opts)
%
% DESCRIPTION:
%   Return a failure flag/message for invalid time, ODE failure
%   or nonfinite augmented output instead of advancing the Newton solve.
%
% INPUT:
%   x0_nd     [6 x 1] SE rotating state [x;y;z;vx;vy;vz], normalized.
%   tf_nd     Positive final time [TU_SE]; initial time is zero.
%   mu_se_nd  Secondary/total mass ratio [-].
%   opts      Struct with ode_rel_tol and ode_abs_tol.
%
% OUTPUT:
%   ok        Logical propagation status.
%   data      Struct with xf_half_nd [6 x 1] and Phi_half [6 x 6],
%               mapping initial-state to terminal-state perturbations.
%   msg       Empty text on success, failure label otherwise.

    ok = false;
    msg = '';
    data = struct('xf_half_nd', [], 'Phi_half', []);

    if ~(isfinite(tf_nd) && tf_nd > 0)
        msg = 'invalid_half_period';
        return;
    end

    Phi0 = eye(6);
    x0_ext = [x0_nd(:); Phi0(:)];
    rhs_stm = @(t, x) rhs_SE3BP_3D_STM(t, x, mu_se_nd);
    ode_opts = odeset('RelTol', opts.ode_rel_tol, 'AbsTol', opts.ode_abs_tol);

    try
        [~, xx_ext] = ode113(rhs_stm, [0, tf_nd], x0_ext, ode_opts);
    catch
        msg = 'half_period_propagation_failed';
        return;
    end

    if isempty(xx_ext) || size(xx_ext, 2) ~= 42 || any(~isfinite(xx_ext(:)))
        msg = 'half_period_propagation_invalid_output';
        return;
    end

    xf_half_nd = xx_ext(end, 1:6).';
    Phi_half = reshape(xx_ext(end, 7:end), 6, 6);

    if any(~isfinite(xf_half_nd)) || any(~isfinite(Phi_half(:)))
        msg = 'half_period_non_finite_terminal';
        return;
    end

    ok = true;
    data.xf_half_nd = xf_half_nd;
    data.Phi_half = Phi_half;
end
