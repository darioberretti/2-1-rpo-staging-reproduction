function [dxdt] = rhs_CR3BP_2D(~,xx, mu)
%RHS_CR3BP_2D Evaluate the 2D CR3BP state derivative.
%
% PROTOTYPE:
%   dxdt = rhs_CR3BP_2D(t, xx, mu)
%
% DESCRIPTION:
%   Evaluate normalized dynamics in the rotating barycentric frame, with
%   primaries at x = -mu and x = 1 - mu. Positions use DU, velocities VU,
%   and derivatives are taken with respect to time in TU.
%
% INPUT:
%   t     Scalar time [TU]; unused because the CR3BP is autonomous.
%   xx    [4 x 1] State [x; y; vx; vy].
%   mu    Secondary-to-total mass ratio [-].
%
% OUTPUT:
%   dxdt  [4 x 1] State derivative, in input order.

    % Read rotating-frame coordinates and velocities.
    x  = xx(1);
    y  = xx(2);
    vx = xx(3);
    vy = xx(4);

    % Compute distances to the two primaries.
    r1 = sqrt((x + mu)^2 + y^2);
    r2 = sqrt((x + mu - 1)^2 + y^2);

    % Evaluate the effective-potential gradient.
    dUdx = x - (1-mu)/r1^3*(mu+x) + mu/r2^3*(1-mu-x);
    dUdy = y - (1-mu)/r1^3*y - mu/r2^3*y;

    % Assemble the state derivative in the documented component order.
    dxdt = zeros(4,1);

    dxdt(1:2) = xx(3:4);
    dxdt(3)   = dUdx + 2*vy;
    dxdt(4)   = dUdy - 2*vx;
    
end
