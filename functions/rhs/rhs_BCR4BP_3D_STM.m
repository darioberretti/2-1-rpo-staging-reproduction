function [dxdt] = rhs_BCR4BP_3D_STM(t, xx, mu, thetaInitial)
%RHS_BCR4BP_3D_STM Evaluate normalized 3D bicircular state/STM dynamics.
%
% PROTOTYPE:
%   dxdt = rhs_BCR4BP_3D_STM(t, xx, mu, thetaInitial)
%
% DESCRIPTION:
%   Use the rotating Earth-Moon barycentric frame and the frozen solar
%   mass, radius, and angular rate. Solar phase evolves as
%   thetaInitial + ws*t; positions are in DU, velocities in VU.
%   The column-major STM satisfies dPhi/dt = (df/dx) Phi and maps
%   perturbations from the initial epoch to the current state.
%
% INPUT:
%   t             Scalar normalized time [TU].
%   xx            [42 x 1] Augmented state [state; Phi(:)].
%                   Physical state order is [x; y; z; vx; vy; vz].
%   mu            Secondary-to-total Earth-Moon mass ratio [-].
%   thetaInitial  Initial Sun phase in the rotating frame [rad].
%
% OUTPUT:
%   dxdt          [42 x 1] State and column-major STM derivatives
%                   with respect to normalized time in TU.

    x_state = xx(1:6);
    Phi = reshape(xx(7:end), 6, 6);

    [f_state, A] = bcr4bp3d_f_and_A(t, x_state, mu, thetaInitial);
    PhiDot = A * Phi;

    dxdt = zeros(42,1);
    dxdt(1:6) = f_state;
    dxdt(7:end) = PhiDot(:);
end

function [f, A] = bcr4bp3d_f_and_A(t, x_state, mu, thetaInitial)
    %BCR4BP3D_F_AND_A Unpack state.
    x = x_state(1); y = x_state(2); z = x_state(3);
    vx = x_state(4); vy = x_state(5); vz = x_state(6);

    % BCR4BP constants in Earth-Moon normalized units.
    ms = 3.28900541e5;
    rho = 3.88811143e2;
    ws = -9.25195985e-1;

    theta0 = 0;
    if nargin > 3 && ~isempty(thetaInitial)
        theta0 = thetaInitial;
    end
    thetaSun = theta0 + ws * t;

    % Sun position in rotating Earth-Moon frame.
    xs = rho * cos(thetaSun);
    ys = rho * sin(thetaSun);
    zs = 0;

    % Relative vectors to primaries and Sun.
    d1 = [x + mu;     y;      z     ];
    d2 = [x + mu - 1; y;      z     ];
    d3 = [x - xs;     y - ys; z - zs];

    r1 = norm(d1);
    r2 = norm(d2);
    r3 = norm(d3);

    r1_3 = r1^3;
    r2_3 = r2^3;
    r3_3 = r3^3;

    % Effective potential gradient components.
    Ux = x ...
       - (1 - mu) * d1(1) / r1_3 ...
       - mu       * d2(1) / r2_3 ...
       - ms       * d3(1) / r3_3 ...
       - ms / rho^2 * cos(thetaSun);

    Uy = y ...
       - (1 - mu) * d1(2) / r1_3 ...
       - mu       * d2(2) / r2_3 ...
       - ms       * d3(2) / r3_3 ...
       - ms / rho^2 * sin(thetaSun);

    Uz = ...
       - (1 - mu) * d1(3) / r1_3 ...
       - mu       * d2(3) / r2_3 ...
       - ms       * d3(3) / r3_3;

    % Accelerations.
    ax = 2 * vy + Ux;
    ay = -2 * vx + Uy;
    az = Uz;

    % Dynamics vector f(x,t).
    f = [vx; vy; vz; ax; ay; az];

    % Hessian of U: point-mass contributions + centrifugal terms.
    H = zeros(3,3);
    H = H + (1 - mu) * hess_m_over_r(d1, r1);
    H = H + mu       * hess_m_over_r(d2, r2);
    H = H + ms       * hess_m_over_r(d3, r3);

    H(1,1) = H(1,1) + 1;
    H(2,2) = H(2,2) + 1;

    % Assemble A = df/dx.
    A = zeros(6,6);
    A(1:3,4:6) = eye(3);
    A(4:6,1:3) = H;

    % Coriolis partials.
    A(4,5) = 2;
    A(5,4) = -2;
end

function H = hess_m_over_r(d, r)
    %HESS_M_OVER_R Hessian of 1/r.
    I = eye(3);
    r3 = r^3;
    r5 = r^5;
    H = 3 * (d * d.') / r5 - I / r3;
end
