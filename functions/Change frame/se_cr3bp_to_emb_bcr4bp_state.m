function [x_em_nd, info] = se_cr3bp_to_emb_bcr4bp_state(x_se_nd, theta_sun_rad, varargin)
%SE_CR3BP_TO_EMB_BCR4BP_STATE Convert SE-CR3BP synodic state to EM/BCR4BP state.
%
% PROTOTYPE:
%   [x_em_nd, info] = se_cr3bp_to_emb_bcr4bp_state(x_se_nd, theta_sun_rad)
%   [x_em_nd, info] = se_cr3bp_to_emb_bcr4bp_state(..., Name, Value, ...)
%
% DESCRIPTION:
%   Map SE to EM normalized rotating states using coplanar frame geometry,
%   length scales and angular rates. Positions use the corresponding DU;
%   velocities use DU/TU, with TU = 1/omega. N is the state-row count.
%   Run 04 explicitly selects se_secondary_center=emb to preserve its
%   historical origin; the general converter defaults to earth.
%
% Reference:
%   van der Weg & Vasile (2015), CMDA, Sec. 2.4, Eqs. (6), (12), (13), (14).
%
% The conversion is performed between two *normalized rotating* frames:
% - Sun-Earth (SE) CR3BP synodic frame (input state x_se_nd),
% - Earth-Moon barycentric (EMB) synodic frame used by EMB BCR4BP (output x_em_nd).
%
% INPUT:
%   x_se_nd        [N x 6] state(s) in SE normalized rotating frame:
%                  [x_se, y_se, z_se, vx_se, vy_se, vz_se]
%   theta_sun_rad  [1 x 1] or [N x 1] Sun phase in EMB rotating frame [rad].
%                  Used to infer the relative geometry angle alpha if
%                  'alpha_rad' is not provided (see options).
%
% NAME-VALUE OPTIONS
%   'alpha_rad'                    [ ] or scalar/[N x 1]
%       Relative geometry angle alpha between EM and SE primaries [rad].
%       If provided, it overrides theta-based reconstruction.
%   'alpha_plus_theta_sun_const_rad' [scalar, default pi]
%       Assumed constant C in alpha + theta_sun = C (coplanar model).
%       Used only when 'alpha_rad' is not provided.
%   'l_se_km'                      [scalar, default 1.495979e8]
%   'l_em_km'                      [scalar, default 384400]
%   'omega_se_rad_s'               [scalar, default 1.99098670e-7]
%   'omega_em_rad_s'               [scalar, default 2.66531437e-6]
%   'mu_se_nd'                     [scalar, default 1/(1+3.28900541e5)]
%   'mu_em_nd'                     [scalar, default 0.0121506683]
%   'se_secondary_center'          [char/string, default 'earth']
%       Use 'earth' when x=1-mu_se is Earth, and 'emb' when x=1-mu_se is
%       the Earth-Moon barycenter as in rhs_SEBCR4BP_3D.
%
% OUTPUT:
%   x_em_nd        [N x 6] state(s) in EMB normalized rotating frame:
%                  [x_em, y_em, z_em, vx_em, vy_em, vz_em]
%   info           Struct with alpha_rad/theta_sun_rad [N x 1, rad],
%                  physical constants above and dimensionless length/rate
%                  ratios; state ordering is [x,y,z,vx,vy,vz] in both frames.
%
% NOTES
% - This function assumes coplanar SE and EM rotating frames.
% - x/y conversion uses complex notation from Eq. (6).
% - velocity conversion follows Eq. (12), inverted to EM from SE.
% - z and vz conversion follow Eqs. (13)-(14), inverted to EM from SE.

    % Parse and validate function inputs/options.
    p = inputParser;
    p.FunctionName = mfilename;
    addRequired(p, 'x_se_nd', @(x) isnumeric(x) && isreal(x) && ismatrix(x) && size(x,2) == 6);
    addRequired(p, 'theta_sun_rad', @(x) isnumeric(x) && isreal(x) && (isscalar(x) || isvector(x)));
    addParameter(p, 'alpha_rad', [], @(x) isnumeric(x) && isreal(x) && (isempty(x) || isscalar(x) || isvector(x)));
    addParameter(p, 'alpha_plus_theta_sun_const_rad', pi, @(x) isnumeric(x) && isreal(x) && isscalar(x));
    addParameter(p, 'l_se_km', 1.495979e8, @(x) isnumeric(x) && isreal(x) && isscalar(x) && x > 0);
    addParameter(p, 'l_em_km', 384400, @(x) isnumeric(x) && isreal(x) && isscalar(x) && x > 0);
    addParameter(p, 'omega_se_rad_s', 1.99098670e-7, @(x) isnumeric(x) && isreal(x) && isscalar(x) && x > 0);
    addParameter(p, 'omega_em_rad_s', 2.66531437e-6, @(x) isnumeric(x) && isreal(x) && isscalar(x) && x > 0);
    addParameter(p, 'mu_se_nd', 1 / (1 + 3.28900541e5), @(x) isnumeric(x) && isreal(x) && isscalar(x) && x > 0 && x < 0.5);
    addParameter(p, 'mu_em_nd', 0.0121506683, @(x) isnumeric(x) && isreal(x) && isscalar(x) && x >= 0 && x < 0.5);
    addParameter(p, 'se_secondary_center', 'earth', @(x) ischar(x) || isstring(x));
    parse(p, x_se_nd, theta_sun_rad, varargin{:});

    % Extract required/optional values in normalized column form.
    x_se_nd = p.Results.x_se_nd;
    theta_sun_rad = p.Results.theta_sun_rad(:);
    alpha_in = p.Results.alpha_rad;

    % Guard against non-finite numerical inputs.
    if any(~isfinite(x_se_nd(:)))
        error('%s:NonFiniteState', mfilename, 'x_se_nd contains NaN/Inf values.');
    end
    if any(~isfinite(theta_sun_rad(:)))
        error('%s:NonFiniteTheta', mfilename, 'theta_sun_rad contains NaN/Inf values.');
    end

    % Broadcast scalar phase inputs to one value per state row.
    n_rows = size(x_se_nd, 1);
    theta_sun_rad = local_expand_scalar_vector(theta_sun_rad, n_rows, 'theta_sun_rad');

    % Reconstruct alpha from thetaSun unless alpha is provided explicitly.
    if isempty(alpha_in)
        alpha_rad = p.Results.alpha_plus_theta_sun_const_rad - theta_sun_rad;
    else
        alpha_rad = local_expand_scalar_vector(alpha_in(:), n_rows, 'alpha_rad');
    end

    % Build nondimensional scaling ratios between SE and EMB systems.
    l_ratio = p.Results.l_se_km / p.Results.l_em_km;
    omega_ratio = p.Results.omega_se_rad_s / p.Results.omega_em_rad_s;
    mu_se_nd = p.Results.mu_se_nd;
    mu_em_nd = p.Results.mu_em_nd;
    secondary_center = lower(string(p.Results.se_secondary_center));
    if secondary_center == "earth"
        eta_secondary_em = -mu_em_nd;
    elseif secondary_center == "emb"
        eta_secondary_em = 0;
    else
        error('%s:InvalidSecondaryCenter', mfilename, ...
            'se_secondary_center must be ''earth'' or ''emb''.');
    end

    % Map planar position from SE complex form to secondary-relative EMB coordinates.
    eta_se = x_se_nd(:,1) + 1i * x_se_nd(:,2);
    eta_secondary_rel_em = l_ratio .* exp(-1i * alpha_rad) .* (eta_se - (1 - mu_se_nd));

    % Shift the secondary-relative coordinates to the EMB barycentric frame.
    eta_em = eta_secondary_rel_em + eta_secondary_em;

    % Map planar velocity using the secondary-relative vector in the frame-rate term.
    deta_se = x_se_nd(:,4) + 1i * x_se_nd(:,5);
    deta_em = l_ratio .* omega_ratio .* exp(-1i * alpha_rad) .* deta_se ...
        - 1i * (1 - omega_ratio) .* eta_secondary_rel_em;

    % Map out-of-plane position/velocity with scale ratios.
    z_em = l_ratio .* x_se_nd(:,3);
    vz_em = l_ratio .* omega_ratio .* x_se_nd(:,6);

    % Assemble EMB output state as [x y z vx vy vz].
    x_em_nd = [real(eta_em), imag(eta_em), z_em, real(deta_em), imag(deta_em), vz_em];

    % Guard against non-finite outputs from invalid constants/geometry.
    if any(~isfinite(x_em_nd(:)))
        error('%s:NonFiniteOutput', mfilename, ...
            'Converted EMB state contains NaN/Inf; verify geometry/constants.');
    end

    % Return conversion metadata for reproducibility/debugging.
    info = struct();
    info.alpha_rad = alpha_rad;
    info.theta_sun_rad = theta_sun_rad;
    info.constants = struct( ...
        'l_se_km', p.Results.l_se_km, ...
        'l_em_km', p.Results.l_em_km, ...
        'omega_se_rad_s', p.Results.omega_se_rad_s, ...
        'omega_em_rad_s', p.Results.omega_em_rad_s, ...
        'mu_se_nd', mu_se_nd, ...
        'mu_em_nd', mu_em_nd, ...
        'se_secondary_center', char(secondary_center), ...
        'l_ratio_se_to_em', l_ratio, ...
        'omega_ratio_se_to_em', omega_ratio, ...
        'alpha_plus_theta_sun_const_rad', p.Results.alpha_plus_theta_sun_const_rad);
end

function vec_out = local_expand_scalar_vector(vec_in, n_rows, vec_name)
%LOCAL_EXPAND_SCALAR_VECTOR Broadcast a scalar or require exactly n_rows entries; return a column.
% vec_name labels size errors; row order follows the input state batch.

    if isscalar(vec_in)
        vec_out = repmat(vec_in, n_rows, 1);
    else
        if numel(vec_in) ~= n_rows
            error('%s:SizeMismatch', mfilename, ...
                '%s must be scalar or have one value per state row (%d).', vec_name, n_rows);
        end
        vec_out = vec_in(:);
    end
end
