function [f_state, A] = rhs_SE3BP_3D_fA(xx, mu_se)
%RHS_SE3BP_3D_FA Evaluate normalized SE CR3BP dynamics and the analytic state Jacobian.
%
% PROTOTYPE:
%   f_state = rhs_SE3BP_3D_fA(xx, mu_se)
%   [f_state, A] = rhs_SE3BP_3D_fA(xx, mu_se)
%
% DESCRIPTION:
%   Use a rotating barycentric frame with Sun at x=-mu_se
%   and secondary at x=1-mu_se. Positions use DU_SE, velocity VU_SE
%   and time TU_SE. Reject ranges <=1e-12 DU_SE from either primary.
%   Dynamics/variations follow Heritier and Howell (AIAA 2012-4666,
%   DOI: 10.2514/6.2012-4666) and Howell (1984, Celestial Mechanics
%   32, 53-71, DOI: 10.1007/BF01358403).
%
% INPUT:
%   xx     [6 elements] State [x,y,z,vx,vy,vz], normalized.
%   mu_se  Scalar secondary/total mass ratio [-].
%
% OUTPUT:
%   f_state  [6 x 1] State derivative with respect to TU_SE.
%   A        [6 x 6] df_state/dxx in the same component order;
%              computed only when requested. The variational equation
%              is Phi_dot=A*Phi.

    x = xx(1);
    y = xx(2);
    z = xx(3);
    vx = xx(4);
    vy = xx(5);
    vz = xx(6);

    % Evaluate primary-relative geometry and reject singular states.
    r1_vec = [x + mu_se; y; z];
    r2_vec = [x - 1 + mu_se; y; z];

    r1 = norm(r1_vec);
    r2 = norm(r2_vec);

    R_MIN = 1e-12;
    if r1 <= R_MIN || r2 <= R_MIN
        error('rhs_SE3BP_3D_fA:Singularity', ...
            'State too close to a primary (r1=%.3e, r2=%.3e).', r1, r2);
    end

    inv_r1_3 = 1 / r1^3;
    inv_r2_3 = 1 / r2^3;

    % Combine effective-potential gradients with rotating-frame Coriolis terms.
    Ux = x - (1 - mu_se) * r1_vec(1) * inv_r1_3 - mu_se * r2_vec(1) * inv_r2_3;
    Uy = y - (1 - mu_se) * r1_vec(2) * inv_r1_3 - mu_se * r2_vec(2) * inv_r2_3;
    Uz =   - (1 - mu_se) * r1_vec(3) * inv_r1_3 - mu_se * r2_vec(3) * inv_r2_3;

    ax = Ux + 2 * vy;
    ay = Uy - 2 * vx;
    az = Uz;

    f_state = [vx; vy; vz; ax; ay; az];

    if nargout < 2
        return;
    end

    % Form the potential Hessian and assemble the optional state Jacobian.
    I3 = eye(3);
    inv_r1_5 = 1 / r1^5;
    inv_r2_5 = 1 / r2^5;

    hess_r1 = 3 * (r1_vec * r1_vec.') * inv_r1_5 - I3 * inv_r1_3;
    hess_r2 = 3 * (r2_vec * r2_vec.') * inv_r2_5 - I3 * inv_r2_3;

    U_hess = (1 - mu_se) * hess_r1 + mu_se * hess_r2;
    U_hess(1,1) = U_hess(1,1) + 1;
    U_hess(2,2) = U_hess(2,2) + 1;

    A = zeros(6, 6);
    A(1:3, 4:6) = I3;
    A(4:6, 1:3) = U_hess;
    A(4,5) = 2;
    A(5,4) = -2;
end
