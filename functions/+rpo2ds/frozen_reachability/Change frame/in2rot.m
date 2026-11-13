function xxRot = in2rot(xx, t, mu)
%IN2ROT Transform Earth-centered inertial states to rotating barycentric axes.
%
% PROTOTYPE:
%   xxRot = in2rot(xx, t, mu)
%
% DESCRIPTION:
%   Undo the rotation through epoch t and the rotating-frame velocity terms,
%   then subtract the Earth-to-barycenter x offset mu. No physical scaling
%   is applied; this is the inverse of rot2in.
%
% INPUT:
%   xx  [N x 4] Rows [xI, yI, vxI, vyI] or [N x 6] rows
%        [xI, yI, zI, vxI, vyI, vzI] in Earth-centered inertial DU/VU.
%   t   Scalar epoch or compatible time column [TU], also frame angle [rad].
%   mu  Secondary-to-total mass ratio [-].
%
% OUTPUT:
%   xxRot  Same-size rotating barycentric states in DU/VU, ordered
%            [x, y, vx, vy] or [x, y, z, vx, vy, vz].

    xI = xx(:,1);
    yI = xx(:,2);

    if size(xx,2) == 4
        vxI = xx(:,3);
        vyI = xx(:,4);
    elseif size(xx,2) == 6
        vxI = xx(:,4);
        vyI = xx(:,5);
    else
        error("The state vectors matrices must be Nx4 or Nx6");
    end

    % --- Time handling (allow scalar or vector) ---
    if isscalar(t), t = t .* ones(size(xI)); end
    c = cos(t);
    s = sin(t);

    % --- Positions: [x+mu; y] = R(-t)*[xI; yI] ---
    x_plus_mu =  xI .* c + yI .* s;
    y_rot     = -xI .* s + yI .* c;
    x_rot     = x_plus_mu - mu;

    % --- Velocities: [vx - y; vy + x + mu] = R(-t)*[vxI; vyI] ---
    a =  vxI .* c + vyI .* s;         % = vx - y
    b = -vxI .* s + vyI .* c;         % = vy + x + mu

    vx_rot = a + y_rot;
    vy_rot = b - x_rot - mu;

    % --- Assemble output ---
    if size(xx,2) == 4
        xxRot = [x_rot, y_rot, vx_rot, vy_rot];
    else
        % z, vz unchanged (rotation is in the xy-plane)
        xxRot = [x_rot, y_rot, xx(:,3), vx_rot, vy_rot, xx(:,6)];
    end
end
