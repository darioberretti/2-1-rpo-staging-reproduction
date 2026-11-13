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

    % Read rotating-frame coordinates and velocities.
    x  = xx(1);
    y  = xx(2);
    vx = xx(3);
    vy = xx(4);

    % Recover the column-major STM at the current epoch.
    Phi = reshape(xx(5:end),4,4);

    % Compute distances to the two primaries.
    r1 = sqrt((x + mu)^2 + y^2);
    r2 = sqrt((x + mu - 1)^2 + y^2);

    % Evaluate the effective-potential gradient.
    dUdx = x - (1-mu)/r1^3*(mu+x) + mu/r2^3*(1-mu-x);
    dUdy = y - (1-mu)/r1^3*y - mu/r2^3*y;

    % Evaluate the planar dynamics Jacobian df/dx.
    dfdx = [0, 0, 1, 0;
            0, 0, 0, 1;
            (mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) + (3*mu*(2*mu + 2*x - 2)*(mu + x - 1))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*(2*mu + 2*x)*(mu + x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2)) + 1,   (3*mu*y*(mu + x - 1))/((mu + x - 1)^2 + y^2)^(5/2) - (3*y*(mu + x)*(mu - 1))/((mu + x)^2 + y^2)^(5/2),  0, 2;
            (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2)), (mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2)^(5/2) + 1, -2, 0];                                                                                                    
                                                                                                                    
    % Propagate perturbations through dPhi/dt = (df/dx) Phi.
    Phidot = dfdx*Phi;

    % Assemble the state derivative in the documented component order.
    dxdt = zeros(20,1);

    dxdt(1:2) = xx(3:4);
    dxdt(3)   = dUdx + 2*vy;
    dxdt(4)   = dUdy - 2*vx;
    dxdt(5:end) = Phidot(:);
    
end
