function [dxdt] = rhs_BCR4BP_3D(t, xx, mu, thetaInitial)
%RHS_BCR4BP_3D Evaluate normalized 3D bicircular state dynamics.
%
% PROTOTYPE:
%   dxdt = rhs_BCR4BP_3D(t, xx, mu, thetaInitial)
%
% DESCRIPTION:
%   Use the rotating Earth-Moon barycentric frame and the frozen solar
%   mass, radius, and angular rate. Solar phase evolves as
%   thetaInitial + ws*t; positions are in DU, velocities in VU.
%
% INPUT:
%   t             Scalar normalized time [TU].
%   xx            [6 x 1] State [x; y; z; vx; vy; vz].
%                   Physical state order is [x; y; z; vx; vy; vz].
%   mu            Secondary-to-total Earth-Moon mass ratio [-].
%   thetaInitial  Initial Sun phase in the rotating frame [rad].
%
% OUTPUT:
%   dxdt          [6 x 1] State derivative
%                   with respect to normalized time in TU.

    x  = xx(1);
    y  = xx(2);
    z  = xx(3);
    vx = xx(4);
    vy = xx(5);

    % BCR4BP constants in EM normalized units
    ms  = 3.28900541e5;   % scaled mass of Sun
    rho = 3.88811143e2;   % Sun-(Earth+Moon barycenter) distance
    ws  = -9.25195985e-1; % Sun angular velocity in EM rotating frame

    thetaSun = thetaInitial + ws * t;
    cTh = cos(thetaSun);
    sTh = sin(thetaSun);

    % Distances
    r1 = sqrt((x + mu)^2     + y^2 + z^2);
    r2 = sqrt((x + mu - 1)^2 + y^2 + z^2);
    r3 = sqrt((x - rho*cTh)^2 + (y - rho*sTh)^2 + z^2);

    % Potential derivatives
    dUdx = x ...
        - (1 - mu) * (x + mu) / r1^3 ...
        + mu * (1 - mu - x) / r2^3 ...
        - ms * (x - rho*cTh) / r3^3 ...
        - ms * cTh / rho^2;

    dUdy = y ...
        - (1 - mu) * y / r1^3 ...
        - mu * y / r2^3 ...
        - ms * (y - rho*sTh) / r3^3 ...
        - ms * sTh / rho^2;

    dUdz = - (1 - mu) * z / r1^3 ...
           - mu * z / r2^3 ...
           - ms * z / r3^3;

    % Dynamics
    dxdt = zeros(6,1);
    dxdt(1:3) = xx(4:6);
    dxdt(4) = dUdx + 2 * vy;
    dxdt(5) = dUdy - 2 * vx;
    dxdt(6) = dUdz;
end

