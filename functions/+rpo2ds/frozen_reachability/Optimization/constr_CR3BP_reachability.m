function [c, ceq, Jc, Jceq] = constr_CR3BP_reachability(yy, constants)
%CONSTR_CR3BP_REACHABILITY Evaluate the frozen fixed-phase reachability constraints.
%
% PROTOTYPE:
%   [c, ceq, Jc, Jceq] = constr_CR3BP_reachability(yy, constants)
%
% DESCRIPTION:
%   Enforce multiple-shooting continuity, departure impulse/state, primary
%   periapsis conditions, outgoing v-infinity direction, time limits, and
%   total normalized maneuver budget. Use the historical node-time layout
%   and derivative implementation without substituting the shared core.
%
% INPUT:
%   yy         [10*N + 1 x 1] Node blocks [X(6); U(3)], N absolute times
%              [TU], and final v-infinity magnitude [km/s]. X uses DU/VU.
%   constants  Legacy vector: mu, rminEarth, rminMoon, muE, target direction
%              (3), departure state (6), DU, VU, N, maxTOF, Moon-periapsis
%              index, dv_budget, then optional case-2b node indices.
%              Radii/time/budget use DU/TU/VU; muE uses km^3/s^2,
%              DU uses km and VU uses km/s.
%
% OUTPUT:
%   c, ceq     Inequality (<= 0) and equality (= 0) columns.
%   Jc, Jceq   Optional dense [numel(yy) x number of constraints] Jacobians.

mu = constants(1);
rminEarth = constants(2);
rminMoon = constants(3);
muE = constants(4);
vinfE_x = constants(5); % Direction of versor
vinfE_y = constants(6); % Direction of versor
vinfE_z = constants(7); % Direction of versor
x_ini = constants(8);
y_ini = constants(9);
z_ini = constants(10);
vx_ini = constants(11);
vy_ini = constants(12);
vz_ini = constants(13);
DU = constants(14);
VU = constants(15);
n = constants(16);
maxTOF = constants(17);
idx_peri_moon = constants(18);
dv_budget = constants(19);

if length(constants) > 19
    I_2b = constants(19+1:end);   
else
    I_2b = [];
end

% Output initialization (TO BE TUNED)
ceq = zeros(6*n+5, 1);  
c = zeros(n+5, 1);

if nargout > 2
    Jceq = zeros(size(ceq,1), 10*n+1);
    Jc = zeros(size(c,1), 10*n+1);
end

% Inequality condition on time interval
t = yy( 9*n+1 : 10*n ) ;
c(1:n-1) = t(1:n-1) - t(2:n) ;  % t(i+1) > t(i)
c(n+2) = t(end) - t(1) - maxTOF; % maximum TOF

if nargout > 2
    % derivatives of time interval
    for i = 1 : n-1
        Jc(i, 9*n+i) = 1 ;
        Jc(i, 9*n+i+1) = -1 ;
    end
    Jc(n+2,9*n+1) = -1;
    Jc(n+2,10*n) = 1;
end

% Inequality condition on deltav magnitude
deltavMag = 0;
for j = 1:n
    base = 9*(j-1);
    u = yy(base + 7);
    w = yy(base + 8);
    s = yy(base + 9);
    deltavMag = deltavMag + (u^2 + w^2 + s^2);
    if nargout > 2
        Jc(n+5, base + 7) = 2*u;
        Jc(n+5, base + 8) = 2*w;
        Jc(n+5, base + 9) = 2*s;
    end
end

c(n+5) = deltavMag - dv_budget;

% initial boundary conditions and their Jacobian

% Extract initial state variables
x1 = yy(1) ; y1 = yy(2) ; z1 = yy(3) ;
vx1 = yy(4) ; vy1 = yy(5) ; vz1 = yy(6) ;
u1 = yy(7) ; w1 = yy(8) ; s1 = yy(9) ;

% Initial equality constraints (from 6*(n-1) to 6*n)
ceq(6*n-5) = x1 - x_ini ;
ceq(6*n-4) = y1 - y_ini ;
ceq(6*n-3) = z1 - z_ini ;

% Check for special 2b case
check1 = length(intersect(I_2b,1)) ;

% Connect velocities through regularized variables (add impulse)
if ~check1  
    ceq(6*n-2) = (vx1 - vx_ini) - (u1^2 - w1^2 - s1^2) ;
    ceq(6*n-1) = (vy1 - vy_ini) - (2 * u1 * w1) ;
    ceq(6*n) = (vz1 - vz_ini) - (2 * u1 * s1) ;
else  
    ceq(6*n-2) = (vx1 - vx_ini) - (2 * u1 * w1) ;
    ceq(6*n-1) = (vy1 - vy_ini) - (u1^2 - w1^2 - s1^2) ;
    ceq(6*n) = (vz1 - vz_ini) - (2 * u1 * s1) ;
end

% Compute jacobians (OK)
if nargout > 2 % if necessary
    
    Jceq(6*n-5, 1) = 1 ;
    Jceq(6*n-4, 2) = 1 ;
    Jceq(6*n-3, 3) = 1 ;
    
    if ~check1  
        Jceq(6*n-2, 4) = 1 ;
        Jceq(6*n-2, 7) = -2 * u1 ;
        Jceq(6*n-2, 8) = 2 * w1 ;
        Jceq(6*n-2, 9) = 2 * s1 ;
        
        Jceq(6*n-1, 5) = 1 ;
        Jceq(6*n-1, 7) = -2 * w1 ;
        Jceq(6*n-1, 8) = -2 * u1 ;
        
        Jceq(6*n, 6) = 1 ;
        Jceq(6*n, 7) = -2 * s1 ;
        Jceq(6*n, 9) = -2 * u1 ;
    else  
        Jceq(6*n-2, 4) = 1 ;
        Jceq(6*n-2, 7) = -2 * w1 ;
        Jceq(6*n-2, 8) = -2 * u1 ;
        
        Jceq(6*n-1, 5) = 1 ;
        Jceq(6*n-1, 7) = -2 * u1 ;
        Jceq(6*n-1, 8) = 2 * w1 ;
        Jceq(6*n-1, 9) = 2 * s1 ;
        
        Jceq(6*n, 6) = 1 ;
        Jceq(6*n, 7) = -2 * s1 ;
        Jceq(6*n, 9) = -2 * u1 ;
    end
end

% Constraints for the intermediate arcs

for i = 1 : n-1

    % Check for special 2b case
    check = length(intersect(I_2b,i+1)) ;  
    
    xi = yy(9*(i-1)+1) ; % current state (x)
    yi = yy(9*(i-1)+2) ; % current state (y)
    zi = yy(9*(i-1)+3) ; % current state (z)
    vxi = yy(9*(i-1)+4) ; % current state (vx)
    vyi = yy(9*(i-1)+5) ; % current state (vy)
    vzi = yy(9*(i-1)+6) ; % current state (vz)
    Xi = [xi, yi, zi, vxi, vyi, vzi]' ;
    
    xip1 = yy(9*i+1) ; % next state (x)
    yip1 = yy(9*i+2) ; % next state (y)
    zip1 = yy(9*i+3) ; % next state (z)
    vxip1 = yy(9*i+4) ; % next state (vx)
    vyip1 = yy(9*i+5) ; % next state (vy)
    vzip1 = yy(9*i+6) ; % next state (vz)
    uip1 = yy(9*i+7) ; % next state (u)
    wip1 = yy(9*i+8) ; % next state (w)
    sip1 = yy(9*i+9) ; % next state (s)

    ti = t(i) ;
    tip1 = t(i+1) ;

    % Compute constraints
    if nargout < 3  % not computing Jacobian (STM not needed)

        % Propagation from a node to the next
        dimensions = '3D';
        framework = 'CR3BP';
        stm_computation = 'false'; % not computing STM
        x_temp  = propagate(Xi,ti,tip1,mu,dimensions,framework,stm_computation);
        
        % Equality constraints (continuity on position and velocity)
        % (from 1 to 3*(n-1))
        ceq(3*(i-1)+1) = x_temp(1) - xip1 ;
        ceq(3*(i-1)+2) = x_temp(2) - yip1 ;
        ceq(3*(i-1)+3) = x_temp(3) - zip1 ;

        % Connect velocities through regularized variables (add impulse)
        % (from 3*(n-1)+1 to 6*(n-1))
        if ~check  
            ceq(3*n-5+3*i) = (vxip1 - x_temp(4)) - (uip1^2 - wip1^2 - sip1^2) ;
            ceq(3*n-4+3*i) = (vyip1 - x_temp(5)) - (2 * uip1 * wip1) ;
            ceq(3*n-3+3*i) = (vzip1 - x_temp(6)) - (2 * uip1 * sip1) ;
        else  
            ceq(3*n-5+3*i) = (vxip1 - x_temp(4)) - (2 * uip1 * wip1) ;
            ceq(3*n-4+3*i) = (vyip1 - x_temp(5)) - (uip1^2 - wip1^2 - sip1^2) ;
            ceq(3*n-3+3*i) = (vzip1 - x_temp(6)) - (2 * uip1 * sip1) ;
        end

    else % computing Jacobian (STM needed)

        % Propagation from a node to the next
        dimensions = '3D';
        framework = 'CR3BP';
        stm_computation = 'true'; % computing STM
        [x_temp, PHI] = propagate(Xi,ti,tip1,mu,dimensions,framework,stm_computation);
        
        % position defect
        ceq(3*(i-1)+1) = x_temp(1) - xip1 ;
        ceq(3*(i-1)+2) = x_temp(2) - yip1 ;
        ceq(3*(i-1)+3) = x_temp(3) - zip1 ;

        % relationship between velocity and regularized variables
        if ~check  
            ceq(3*n-5+3*i) = (vxip1 - x_temp(4)) - (uip1^2 - wip1^2 - sip1^2) ;
            ceq(3*n-4+3*i) = (vyip1 - x_temp(5)) - (2 * uip1 * wip1) ;
            ceq(3*n-3+3*i) = (vzip1 - x_temp(6)) - (2 * uip1 * sip1) ;
        else 
            ceq(3*n-5+3*i) = (vxip1 - x_temp(4)) - (2 * uip1 * wip1) ;
            ceq(3*n-4+3*i) = (vyip1 - x_temp(5)) - (uip1^2 - wip1^2 - sip1^2) ;
            ceq(3*n-3+3*i) = (vzip1 - x_temp(6)) - (2 * uip1 * sip1) ;
        end

        % Values at i-th mesh
        fi = rhs_CR3BP_3D(ti,Xi, mu); % ti not needed

        % Values after integrating i-th mesh
        f_temp = rhs_CR3BP_3D(tip1,x_temp, mu); % ti not needed

        % position defect (OK)
        Jceq(3*(i-1)+1: 3*i, 9*(i-1)+1: 9*(i-1)+6) = PHI(1:3, 1:6) ;   % 位置だけのdefectでも，速度に依存することに注意！！！
        Jceq(3*(i-1)+1: 3*i, 9*i+1: 9*i+3) = - eye(3) ;

        % defect derivatives wrt t (OK)
        Jceq(3*(i-1)+1:3*i, 9*n+i) = - PHI(1:3, 1:6) * fi(1:6, 1) ;  % cc_eq derivatives wrt tt(i)
        Jceq(3*(i-1)+1:3*i, 9*n+i+1) = f_temp(1:3, 1) ;  % cc_eq derivatives wrt tt(i+1)

        % relationship between velocity and regularized variables (OK)
        Jceq(3*n-5+3*i:3*n-3+3*i, 9*i+4:9*i+6) = eye(3) ;
        Jceq(3*n-5+3*i:3*n-3+3*i, 9*(i-1)+1:9*(i-1)+6) = -PHI(4:6, 1:6) ;
        Jceq(3*n-5+3*i:3*n-3+3*i, 9*n+i) = PHI(4:6, 1:6) * fi(1:6, 1) ;
        Jceq(3*n-5+3*i:3*n-3+3*i, 9*n+i+1) = -f_temp(4:6, 1) ;
        
        if check == 0 % (rows from 3*n-2 to 6*n-6) (OK)
            Jceq(3*n-5+3*i:3*n-3+3*i, 9*i+7:9*i+9) = 2 * [-uip1, wip1, sip1; -wip1, -uip1, 0; -sip1, 0, -uip1] ;  
        else 
            Jceq(3*n-5+3*i:3*n-3+3*i, 9*i+7:9*i+9) = 2 * [-wip1, -uip1, 0; -uip1, wip1, sip1; -sip1, 0, -uip1] ;  
        end
        
    end  

end

% Moon flyby - minimum distance of rp from Moon

X_fb = yy(9*idx_peri_moon-8:9*idx_peri_moon-3)' ;

% r1, r1dot, r1ddot calculation
[r2,r2dot,r2ddot] = get_r2ddot(X_fb, mu) ;

% Constraints to ensure the propagation ends at the perigee of the flyby
c(n+3) = rminMoon^2 - r2^2 ;  
ceq(6*n+5) = r2dot ;  % r2dot = 0 
c(n+4) = - r2ddot ;  % r2ddot > 0 

if nargout > 2
    [dr2dot_dt, dr2dot_dx, dr2dot_dy, dr2dot_dz, ...
          dr2dot_dvx, dr2dot_dvy, dr2dot_dvz, ...
          dr2ddot_dt, dr2ddot_dx, dr2ddot_dy, dr2ddot_dz, ...
          dr2ddot_dvx, dr2ddot_dvy, dr2ddot_dvz] = get_dr2ddot_dtX(X_fb, constants);

    Jc(n+3,9*idx_peri_moon-8) = -2*(X_fb(1)-1+mu) ;
    Jc(n+3,9*idx_peri_moon-7) = -2*X_fb(2) ;
    Jc(n+3,9*idx_peri_moon-6) = -2*X_fb(3) ;
    Jceq(6*n+5,9*n+idx_peri_moon) = dr2dot_dt ;
    Jceq(6*n+5,9*idx_peri_moon-8) = dr2dot_dx ;
    Jceq(6*n+5,9*idx_peri_moon-7) = dr2dot_dy ;
    Jceq(6*n+5,9*idx_peri_moon-6) = dr2dot_dz ;
    Jceq(6*n+5,9*idx_peri_moon-5) = dr2dot_dvx ;
    Jceq(6*n+5,9*idx_peri_moon-4) = dr2dot_dvy ;
    Jceq(6*n+5,9*idx_peri_moon-3) = dr2dot_dvz ;
    
    Jc(n+4,9*n+idx_peri_moon) = -dr2ddot_dt ;
    Jc(n+4,9*idx_peri_moon-8) = -dr2ddot_dx ; % PROBLEM
    Jc(n+4,9*idx_peri_moon-7) = -dr2ddot_dy ; % PROBLEM
    Jc(n+4,9*idx_peri_moon-6) = -dr2ddot_dz ; % PROBLEM
    Jc(n+4,9*idx_peri_moon-5) = -dr2ddot_dvx ;
    Jc(n+4,9*idx_peri_moon-4) = -dr2ddot_dvy ;
    Jc(n+4,9*idx_peri_moon-3) = -dr2ddot_dvz ;
end
% Final boundary conditions and their Jacobian (ceq from 6*n+1 to 6*n+4)

Xn = yy(9*n-8:9*n-3)' ;
tn = yy(10*n) ;
vInfMag = yy(10*n+1);

% r1, r1dot, r1ddot calculation
[r1,r1dot,r1ddot] = get_r1ddot(Xn, mu) ;

% Constraints to ensure the propagation ends at the perigee of the flyby
c(n+1) = rminEarth^2 - r1^2 ;  
ceq(6*n+1) = r1dot ;  % r1dot = 0 
c(n) = - r1ddot ;  % r1ddot > 0 

% Compute vinf_out
Xn_in = rot2in(Xn,tn,mu);
Xn_in(1:3) = Xn_in(1:3) * DU;
Xn_in(4:6) = Xn_in(4:6) * VU;
[~, ~, vinf_out_DVaft] = earth_flyby_for_optimization(Xn_in(1:3), Xn_in(4:6), muE) ;

%　constraint for vinf_out
ceq(6*n+2) = vinf_out_DVaft(1) - vinfE_x * vInfMag;
ceq(6*n+3) = vinf_out_DVaft(2) - vinfE_y * vInfMag;
ceq(6*n+4) = vinf_out_DVaft(3) - vinfE_z * vInfMag;

% compute Jacobian
if nargout > 2
    
    [dr1dot_dt,dr1dot_dx,dr1dot_dy,dr1dot_dz,dr1dot_dvx,dr1dot_dvy,dr1dot_dvz, ...,
        dr1ddot_dt,dr1ddot_dx,dr1ddot_dy,dr1ddot_dz,dr1ddot_dvx,dr1ddot_dvy,dr1ddot_dvz] = get_dr1ddot_dtX(Xn,constants) ;
    
    Jc(n+1,9*n-8) = -2*(Xn(1)+mu) ;
    Jc(n+1,9*n-7) = -2*Xn(2) ;
    Jc(n+1,9*n-6) = -2*Xn(3) ;
    Jceq(6*n+1,10*n) = dr1dot_dt ;
    Jceq(6*n+1,9*n-8) = dr1dot_dx ;
    Jceq(6*n+1,9*n-7) = dr1dot_dy ;
    Jceq(6*n+1,9*n-6) = dr1dot_dz ;
    Jceq(6*n+1,9*n-5) = dr1dot_dvx ;
    Jceq(6*n+1,9*n-4) = dr1dot_dvy ;
    Jceq(6*n+1,9*n-3) = dr1dot_dvz ;
    
    Jc(n,10*n) = -dr1ddot_dt ;
    Jc(n,9*n-8) = -dr1ddot_dx ; % PROBLEM
    Jc(n,9*n-7) = -dr1ddot_dy ; % PROBLEM
    Jc(n,9*n-6) = -dr1ddot_dz ; % PROBLEM
    Jc(n,9*n-5) = -dr1ddot_dvx ;
    Jc(n,9*n-4) = -dr1ddot_dvy ;
    Jc(n,9*n-3) = -dr1ddot_dvz ;
    
    % use central difference
    dvinfdXn_block = dvinfdXn (Xn, tn, constants) ;
    Jceq(6*n+2:6*n+4, 9*n-8:9*n-3)  = dvinfdXn_block(:,1:6) ; % Some issues here but not too bad
    Jceq(6*n+2:6*n+4, 10*n) = dvinfdXn_block(:,7);
    Jceq(6*n+2, 10*n+1) = -vinfE_x;
    Jceq(6*n+3, 10*n+1) = -vinfE_y;
    Jceq(6*n+4, 10*n+1) = -vinfE_z;
    
end

if nargout > 2

    Jceq = Jceq' ;
    Jc = Jc' ;

end
end
%
function dvinfdXn_block = dvinfdXn(Xn, tn, constants)
%DVINFDXN
% Computes the Jacobian of the outgoing hyperbolic excess velocity vector
% v_inf_out with respect to the final after-impulse state Xn and time tn.
%
%   dvinfdXn_block(i,j) = ∂(v_inf_out_i) / ∂(Xn_j)
%
% using central finite differences.
%
% INPUT:
%   Xn        : 6×1 state vector at the final node (rotating CR3BP frame)
%               Xn = [x; y; z; vx; vy; vz]
%
%   tn        : scalar time at the final node
%
%   constants : parameter vector
%               constants(1) = mu  (CR3BP mass parameter)
%
% OUTPUT:
%   dvinfdXn_block : 3×6 matrix of partial derivatives

mu = constants(1);
muE = constants(4);
DU = constants(14);
VU = constants(15);

% Preallocate Jacobian matrix
dvinfdXn_block = zeros(3,7);

% Optimal finite difference step
delta = sqrt(eps);

for i = 1:6
    % 1) Build perturbed states (Xn + delta e_i) and (Xn - delta e_i)
    Xn_plus  = Xn;
    Xn_minus = Xn;
    Xn_plus(i)  = Xn_plus(i)  + delta;
    Xn_minus(i) = Xn_minus(i) - delta;
    % 2) Convert perturbed states to inertial ECI frame
    Xin_plus  = rot2in(Xn_plus,  tn, mu);
    Xin_plus(1:3) = Xin_plus(1:3) * DU;
    Xin_plus(4:6) = Xin_plus(4:6) * VU;

    Xin_minus = rot2in(Xn_minus, tn, mu);
    Xin_minus(1:3) = Xin_minus(1:3) * DU;
    Xin_minus(4:6) = Xin_minus(4:6) * VU;
    % 3) Compute outgoing v_inf for both perturbed states
    [~, ~, vinf_plus]  = earth_flyby_for_optimization(Xin_plus(1:3),  Xin_plus(4:6),  muE);
    [~, ~, vinf_minus] = earth_flyby_for_optimization(Xin_minus(1:3), Xin_minus(4:6), muE);
    % 4) Central finite difference derivative
    dvinfdXn_block(:,i) = (vinf_plus(:) - vinf_minus(:)) / (2 * delta);

end

% TIME GRADIENTS
% 1) Build perturbed times
tn_plus  = tn + delta;
tn_minus = tn - delta;
% 2) Convert perturbed states to inertial ECI frame
Xin_plus  = rot2in(Xn,  tn_plus, mu);
Xin_plus(1:3) = Xin_plus(1:3) * DU;
Xin_plus(4:6) = Xin_plus(4:6) * VU;

Xin_minus = rot2in(Xn, tn_minus, mu);
Xin_minus(1:3) = Xin_minus(1:3) * DU;
Xin_minus(4:6) = Xin_minus(4:6) * VU;
% 3) Compute outgoing v_inf for both perturbed states
[~, ~, vinf_plus]  = earth_flyby_for_optimization(Xin_plus(1:3),  Xin_plus(4:6),  muE);
[~, ~, vinf_minus] = earth_flyby_for_optimization(Xin_minus(1:3), Xin_minus(4:6), muE);
% 4) Central finite difference derivative
dvinfdXn_block(:,7) = (vinf_plus(:) - vinf_minus(:)) / (2 * delta);

end
