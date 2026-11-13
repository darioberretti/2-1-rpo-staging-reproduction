function [dxdt] = rhs_BCR4BP_2D_STM(t,xx, mu, thetaInitial)
%RHS_BCR4BP_2D_STM Evaluate normalized 2D bicircular state/STM dynamics.
%
% PROTOTYPE:
%   dxdt = rhs_BCR4BP_2D_STM(t, xx, mu, thetaInitial)
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
%   xx            [20 x 1] Augmented state [state; Phi(:)].
%                   Physical state order is [x; y; vx; vy].
%   mu            Secondary-to-total Earth-Moon mass ratio [-].
%   thetaInitial  Initial Sun phase in the rotating frame [rad].
%
% OUTPUT:
%   dxdt          [20 x 1] State and column-major STM derivatives
%                   with respect to normalized time in TU.

    if nargin < 4 || isempty(thetaInitial)
        thetaInitial = 0;
    end

    % Extract variables
    x  = xx(1);
    y  = xx(2);
    vx = xx(3);
    vy = xx(4);
    
    % Other needed variables
    ms = 3.28900541e5;                          % [-] Scaled mass of the Sun
    rho =  3.88811143e2;                        % [-] Scaled Sun-(Earth+Moon) distance
    ws = -9.25195985e-1;                        % [-] Scaled angular velocity of the Sun
    theta0 = thetaInitial;                      % [rad] Initial angle of the sun
    thetaSun = theta0 + ws*t;                   % [rad] Angle of the sun

    % Put PHI in matrix form
    Phi = reshape(xx(5:end),4,4);

    % Compute distances from bodies 1 and 2
    r1 = sqrt((x + mu)^2 + y^2);
    r2 = sqrt((x + mu - 1)^2 + y^2);
    r3 = sqrt((x - rho*cos(thetaSun))^2 + (y - rho*sin(thetaSun))^2);

    % Compute derivative of the potential
    dUdx = x - (1-mu)/r1^3*(mu+x) + mu/r2^3*(1-mu-x) - ms*(x - rho*cos(thetaSun))/r3^3 - ms/rho^2*cos(thetaSun);
    dUdy = y - (1-mu)/r1^3*y - mu/r2^3*y - ms*(y - rho*sin(thetaSun))/r3^3 - ms/rho^2*sin(thetaSun);

    % Assemble the matrix A(t)=dfdx 4x4 matrix
    dfdx = [0, 0,  1, 0;
            0, 0,  0, 1;
           (mu - 1)/((mu + x)^2 + y^2)^(3/2) - ms/((x - rho*cos(thetaSun))^2 + (y - rho*sin(thetaSun))^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) + (3*ms*(x - rho*cos(thetaSun))*(2*x - 2*rho*cos(thetaSun)))/(2*((x - rho*cos(thetaSun))^2 + (y - rho*sin(thetaSun))^2)^(5/2)) + (3*mu*(2*mu + 2*x - 2)*(mu + x - 1))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*(2*mu + 2*x)*(mu + x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2)) + 1,                                                                                                                    (3*ms*(x - rho*cos(thetaSun))*(2*y - 2*rho*sin(thetaSun)))/(2*((x - rho*cos(thetaSun))^2 + (y - rho*sin(thetaSun))^2)^(5/2)) - (3*y*(mu + x)*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*mu*y*(mu + x - 1))/((mu + x - 1)^2 + y^2)^(5/2),  0, 2;
           (3*ms*(y - rho*sin(thetaSun))*(2*x - 2*rho*cos(thetaSun)))/(2*((x - rho*cos(thetaSun))^2 + (y - rho*sin(thetaSun))^2)^(5/2)) + (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2)), (mu - 1)/((mu + x)^2 + y^2)^(3/2) - ms/((x - rho*cos(thetaSun))^2 + (y - rho*sin(thetaSun))^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2)^(5/2) + (3*ms*(y - rho*sin(thetaSun))*(2*y - 2*rho*sin(thetaSun)))/(2*((x - rho*cos(thetaSun))^2 + (y - rho*sin(thetaSun))^2)^(5/2)) + 1, -2, 0];

    % Compute the derivative of the STM
    Phidot = dfdx*Phi;

    % Assemble right-hand side
    dxdt = zeros(20,1);

    dxdt(1:2) = xx(3:4);
    dxdt(3)   = dUdx + 2*vy;
    dxdt(4)   = dUdy - 2*vx;
    dxdt(5:end) = Phidot(:);
    
end
