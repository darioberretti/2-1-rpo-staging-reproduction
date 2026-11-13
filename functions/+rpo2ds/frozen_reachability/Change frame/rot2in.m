function xxInertial = rot2in(xx, t, mu)
%ROT2IN Transform normalized rotating states to Earth-centered inertial axes.
%
% PROTOTYPE:
%   xxInertial = rot2in(xx, t, mu)
%
% DESCRIPTION:
%   Shift the barycentric x coordinate by mu and rotate through t radians.
%   Include the rotating-frame angular-velocity contribution in velocity.
%   This changes frame/origin without applying physical scale factors.
%
% INPUT:
%   xx  [N x 4] Rows [x, y, vx, vy] or [N x 6] [x, y, z, vx, vy, vz].
%        Rotating barycentric positions/velocities use DU/VU.
%   t   Scalar epoch or compatible time column [TU], also the frame angle [rad].
%   mu  Secondary-to-total mass ratio [-].
%
% OUTPUT:
%   xxInertial  Same-size Earth-centered inertial rows in DU/VU.

    x = xx(:,1);
    y = xx(:,2);
    if size(xx,2) == 4
        vx = xx(:,3);
        vy = xx(:,4);
    elseif size(xx,2) == 6
        vx = xx(:,4);
        vy = xx(:,5);
    else
        error("The state vectors matrices must be Nx4 or Nx6");
    end

    % Transform into inertial frame (x+mu --> x if we want RF centered on
    % the EMB: now it is Earth centered)
    xIn = (x + mu) .* cos(t) - y .* sin(t);
    yIn = (x + mu) .* sin(t) + y .* cos(t);
    vxIn = (vx - y) .* cos(t) - (vy + x + mu) .* sin(t);
    vyIn = (vx - y) .* sin(t) + (vy + x +  mu) .* cos(t);

    % Position and velocity vectors
    if size(xx,2) == 4
        xxInertial = [xIn, yIn, vxIn, vyIn];
    elseif size(xx,2) == 6
        xxInertial = [xIn, yIn, xx(:,3), vxIn, vyIn, xx(:,6)];
    end
    
end
