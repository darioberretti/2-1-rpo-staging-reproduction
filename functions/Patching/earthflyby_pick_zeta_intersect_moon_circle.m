function sol = earthflyby_pick_zeta_intersect_moon_circle(vinf_plus, rp, muE, rM_vec, vM_vec, do_warn)
%EARTHFLYBY_PICK_ZETA_INTERSECT_MOON_CIRCLE Choose B-plane azimuths intersecting a circular lunar orbit plane.
%
% PROTOTYPE:
%   sol = earthflyby_pick_zeta_intersect_moon_circle( ...
%       vinf_plus, rp, muE, rM_vec, vM_vec, do_warn)
%
% DESCRIPTION:
%   Find the incoming hyperbola location at the lunar orbital
%   radius, then solve the plane-intersection equation in B-plane
%   azimuth. Return exact solutions when available, a closest azimuth
%   otherwise, or empty geometry when the hyperbola cannot reach the
%   radius. Periapsis states use rp_circle_from_vinfplus.
%
% INPUT:
%   vinf_plus  [3 elements] Outgoing Earth-relative excess velocity [km/s].
%   rp         Positive Earth periapsis radius [km].
%   muE        Earth gravitational parameter [km^3/s^2].
%   rM_vec     [3 elements] Earth-to-Moon position [km], setting orbit radius.
%   vM_vec     [3 elements] Lunar velocity [km/s], setting the orbit normal.
%   do_warn    Optional logical warning flag, default false.
%
% OUTPUT:
%   sol  Struct with has_exact, zeta [1 x K, rad], Rp/Vp [K x 3,
%          Earth-centered periapsis km and km/s], f_in [rad], r_int
%          [K x 3, km], vhat_int/s_minus [K x 3, directions] and
%          plane_resid [1 x K, normal dot incoming radius direction].
%          K is normally 2 for exact solutions or 1 for a closest/
%          degenerate solution; unreachable lunar radius returns K=0.
%          vhat_int retains the historical transverse-direction formula.

    if nargin < 6
        do_warn = false;
    end

    s = vinf_plus(:).';                         % row
    vinf = norm(s);
    s = s / vinf;

    % Hyperbola geometry
    e   = 1 + rp*(vinf^2)/muE;
    if e <= 1, error('Invalid e<=1 (check rp, vinf, muE).'); end
    nu_inf = acos(-1/e);
    vp = sqrt(vinf^2 + 2*muE/rp);               
    h  = rp*vp;
    p  = h^2 / muE;
    rM = norm(rM_vec);
    nMhat = cross(rM_vec, vM_vec);  nMhat = nMhat / norm(nMhat);

    % Feasibility of hitting r=rM on this hyperbola
    cosf = (p/rM - 1)/e;
    if abs(cosf) > 1+1e-12
        sol = struct('has_exact',false,'zeta',[],'Rp',[],'Vp',[], ...
            'f_in',NaN,'r_int',[],'vhat_int',[],'s_minus',[],'plane_resid',[]);
        if do_warn
            warning('No real f at r=rM (|cos f|>1). Hyperbola never reaches rM.');
        end
        return
    end
    cosf = max(-1,min(1,cosf));
    f_in = -acos(cosf);                         % incoming branch
    if ~( -nu_inf < f_in && f_in < 0 )
        % Numerically enforce incoming branch bounds
        f_in = max(-nu_inf+1e-12, min(-1e-12, f_in));
    end
    gamma = f_in - nu_inf;

    % B-plane basis {e1,e2}
    ref = [1 0 0]; if abs(dot(ref,s)) > 0.99, ref = [0 1 0]; end
    e1 = cross(s,ref); e1 = e1 / norm(e1);
    e2 = cross(s,e1);

    % Coefficients for A cos zeta + B sin zeta = K
    A = dot(nMhat, e1);
    B = dot(nMhat, e2);
    c = dot(nMhat, s);
    K = c * cot(gamma);                         % RHS
    R = hypot(A,B);  phi = atan2(B, A);

    solutions = [];
    has_exact = false;

    if R < 1e-14
        % nMhat || s  -> B-plane is perpendicular to nMhat; n·rhat = -c cos(gamma)
        if abs(c*cos(gamma)) < 1e-12
            solutions = phi; has_exact = true;  % any zeta works; pick phi
        else
            % choose arbitrary closest
            solutions = phi; has_exact = false;
        end
    else
        q = K / R;
        if abs(q) <= 1 + 1e-12
            q = max(-1,min(1,q));
            solutions = [phi + acos(q), phi - acos(q)];
            has_exact = true;
        else
            % No exact solution: choose closest (clamp)
            q = max(-1,min(1,q));
            solutions = phi + acos(q);          % single best
            has_exact = false;
        end
    end

    % Build outputs for each solution
    Ksol = numel(solutions);
    Rp = zeros(Ksol,3); Vp = zeros(Ksol,3);
    r_int = zeros(Ksol,3); vhat_int = zeros(Ksol,3);
    s_minus = zeros(Ksol,3); plane_resid = zeros(1,Ksol);

    for k=1:Ksol
        z = solutions(k);
        Bhat = cos(z)*e1 + sin(z)*e2;
        hhat = cross(s, Bhat);  % unit automatically

        % Incoming asymptote direction (Earth flyby)
        delta = 2*asin(1/e);
        s_minus(k,:) =  cos(delta)*s - sin(delta)*Bhat;

        % Intersection direction at r=rM (incoming)
        rhat_fk = -cos(gamma)*s + sin(gamma)*Bhat;
        r_int(k,:) = rM * rhat_fk;
        plane_resid(k) = dot(nMhat, rhat_fk);  % ~0 if exact

        % Velocity direction at intersection (unit, orbital-plane)
        that_p = -sin(nu_inf)*s + cos(nu_inf)*Bhat;
        rhat_p = cross(hhat, that_p);
        that_fk = -sin(f_in)*rhat_p + cos(f_in)*that_p; % unit transverse at f_in
        vhat_int(k,:) = that_fk;

        % Periapsis state via your sampler (Earth-centered)
        [Rp(k,:), Vp(k,:)] = rp_circle_from_vinfplus(vinf_plus, rp, muE, 1, z);
    end

    sol = struct('has_exact',has_exact, 'zeta',solutions, ...
                 'Rp',Rp, 'Vp',Vp, 'f_in',f_in, ...
                 'r_int',r_int, 'vhat_int',vhat_int, ...
                 's_minus',s_minus, 'plane_resid',plane_resid);
end
