function [dxdt] = rhs_CR3BP_2D_STM(~,xx, mu)
%RHS_CR3BP_2D_STM Evaluate the 2D CR3BP state and STM derivatives.
%
% PROTOTYPE:
%   dxdt = rhs_CR3BP_2D_STM(t, xx, mu)
%
% DESCRIPTION:
%   Evaluate normalized dynamics in the rotating barycentric frame, with
%   primaries at x = -mu and x = 1 - mu. Positions use DU, velocities VU,
%   and derivatives are taken with respect to time in TU.
%   The STM satisfies dPhi/dt = (df/dx) Phi and maps perturbations from
%   the initial epoch to the current state. Its 4 x 4 matrix is flattened
%   in MATLAB column-major order.
%
% INPUT:
%   t     Scalar time [TU]; unused because the CR3BP is autonomous.
%   xx    [20 x 1] Augmented state [x; y; vx; vy; Phi(:)].
%   mu    Secondary-to-total mass ratio [-].
%
% OUTPUT:
%   dxdt  [20 x 1] State and column-major STM derivatives, in input order.

    x  = xx(1);
    y  = xx(2);
    vx = xx(3);
    vy = xx(4);

    % Put PHI in matrix form
    Phi = reshape(xx(5:end),4,4);

    % Compute distances from bodies 1 and 2
    r1 = sqrt((x + mu)^2 + y^2);
    r2 = sqrt((x + mu - 1)^2 + y^2);

    % Compute derivative of the potential
    dUdx = x - (1-mu)/r1^3*(mu+x) + mu/r2^3*(1-mu-x);
    dUdy = y - (1-mu)/r1^3*y - mu/r2^3*y;

    % Assemble the matrix A(t)=dfdx 4*4 matrix
    dfdx = [0, 0, 1, 0;
            0, 0, 0, 1;
            (mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) + (3*mu*(2*mu + 2*x - 2)*(mu + x - 1))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*(2*mu + 2*x)*(mu + x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2)) + 1,   (3*mu*y*(mu + x - 1))/((mu + x - 1)^2 + y^2)^(5/2) - (3*y*(mu + x)*(mu - 1))/((mu + x)^2 + y^2)^(5/2),  0, 2;
            (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2)), (mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2)^(5/2) + 1, -2, 0];                                                                                                    
                                                                                                                    
    % Compute the derivative of the STM
    Phidot = dfdx*Phi;

    % Assemble right-hand side
    dxdt = zeros(20,1);

    dxdt(1:2) = xx(3:4);
    dxdt(3)   = dUdx + 2*vy;
    dxdt(4)   = dUdy - 2*vx;
    dxdt(5:end) = Phidot(:);
    
end
