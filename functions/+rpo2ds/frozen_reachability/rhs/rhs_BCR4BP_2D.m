function [dxdt] = rhs_BCR4BP_2D(t,xx, mu, thetaInitial)
%RHS_BCR4BP_2D Evaluate normalized 2D bicircular state dynamics.
%
% PROTOTYPE:
%   dxdt = rhs_BCR4BP_2D(t, xx, mu, thetaInitial)
%
% DESCRIPTION:
%   Use the rotating Earth-Moon barycentric frame and the frozen solar
%   mass, radius, and angular rate. Solar phase evolves as
%   thetaInitial + ws*t; positions are in DU, velocities in VU.
%
% INPUT:
%   t             Scalar normalized time [TU].
%   xx            [4 x 1] State [x; y; vx; vy].
%                   Physical state order is [x; y; vx; vy].
%   mu            Secondary-to-total Earth-Moon mass ratio [-].
%   thetaInitial  Initial Sun phase in the rotating frame [rad].
%
% OUTPUT:
%   dxdt          [4 x 1] State derivative
%                   with respect to normalized time in TU.

    x  = xx(1);
    y  = xx(2);
    vx = xx(3);
    vy = xx(4);
    
    % Other needed variables
    ms = 3.28900541e5;                          % [-] Scaled mass of the Sun
    rho =  3.88811143e2;                        % [-] Scaled Sun-(Earth+Moon) distance
    ws = -9.25195985e-1;                        % [rad/TU] Scaled angular velocity of the Sun
    theta0 = thetaInitial;                      % [rad] Initial angle of the sun
    thetaSun = theta0 + ws*t;                   % [rad] Angle of the sun

    % Compute distances from bodies 1 and 2
    r1 = sqrt((x + mu)^2 + y^2);
    r2 = sqrt((x + mu - 1)^2 + y^2);
    r3 = sqrt((x - rho*cos(thetaSun))^2 + (y - rho*sin(thetaSun))^2);

    % Compute derivative of the potential
    dUdx = x - (1-mu)/r1^3*(mu+x) + mu/r2^3*(1-mu-x) - ms*(x - rho*cos(thetaSun))/r3^3 - ms/rho^2*cos(thetaSun);
    dUdy = y - (1-mu)/r1^3*y - mu/r2^3*y - ms*(y - rho*sin(thetaSun))/r3^3 - ms/rho^2*sin(thetaSun);

    % Assemble right-hand side
    dxdt = zeros(4,1);

    dxdt(1:2) = xx(3:4);
    dxdt(3)   = dUdx + 2*vy;
    dxdt(4)   = dUdy - 2*vx;
    
end
