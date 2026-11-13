function [xf,STMf,STM_hist, xx, tt, tf]  = propagate(x0,t0,tf,mu,dimensions,framework,stm_computation,varargin)
%PROPAGATE Integrate the frozen rotating-frame dynamics and optional STM.
%
% PROTOTYPE:
%   [xf, STMf, STM_hist, xx, tt, tf] = propagate(x0, t0, tf, mu, ...
%       dimensions, framework, stm_computation, evtFlag, evtDir, ...
%       thetaSun0, absTol, relTol)
%
% DESCRIPTION:
%   Use ode113 with the historical positional interface. Automatic STM
%   integration starts from identity and maps initial perturbations to the
%   segment end. This isolated version has no struct/custom-RHS interface.
%
% INPUT:
%   x0               [4 or 6 elements] [x, y, vx, vy] or
%                      [x, y, z, vx, vy, vz] in rotating DU/VU.
%   t0, tf           Initial/requested final epochs [TU].
%   mu               Secondary-to-total mass ratio [-].
%   dimensions       Character vector 2D or 3D.
%   framework        CR3BP, or BCR4BP for the implemented planar branches.
%   stm_computation  Character vector true or false.
%   evtFlag          Optional terminal y = 0 event flag, default false.
%   evtDir           Optional crossing direction -1/0/1, default 0.
%   thetaSun0        Optional initial solar phase [rad], default 0.
%   absTol, relTol   Optional tolerances, both default 1e-12; [] uses defaults.
%
% OUTPUT:
%   xf        [n x 1] Final physical state in DU/VU.
%   STMf      [n x n] Final STM, or [] when disabled.
%   STM_hist  [n x n x Nt] STM history, or [].
%   xx        [Nt x n] State rows with n^2 additional column-major STM
%               columns when enabled. A zero-length interval returns only
%               the physical state row and identity STM.
%   tt, tf    Returned time column and actual final epoch [TU].

    if nargin > 7 && ~isempty(varargin{1})
        evtFlag=varargin{1};
    else
        evtFlag=false;
    end

    if nargin > 8 && ~isempty(varargin{2})
        evtDir=varargin{2};
    else
        evtDir=0;
    end

    if nargin > 9 && ~isempty(varargin{3})
        thetaSun0=varargin{3};
    else
        thetaSun0 = 0;
    end

    if nargin > 10 && ~isempty(varargin{4})
        absTol=varargin{4};
    else
        absTol = 1e-12;
    end

    if nargin > 11 && ~isempty(varargin{5})
        relTol=varargin{5};
    else
        relTol = 1e-12;
    end

    n_states = numel(x0); % number of states
    x0 = x0(:); % verticalize x0

    if strcmpi(stm_computation,'true')
        % Initialize State Transition Matrix at t0
        Phi0 = eye(n_states);
        % Append to initial conditions the conditions for the STM
        x0 = [x0; Phi0(:)];
    end

    switch [dimensions framework stm_computation]
        case ['2D' 'CR3BP' 'true']
            rhs = @(t,x) rhs_CR3BP_2D_STM(t,x,mu);
        case ['2D' 'CR3BP' 'false']
            rhs = @(t,x) rhs_CR3BP_2D(t,x,mu);
        case ['2D' 'BCR4BP' 'true']
            rhs = @(t,x) rhs_BCR4BP_2D_STM(t,x,mu,thetaSun0);
        case ['2D' 'BCR4BP' 'false']
            rhs = @(t,x) rhs_BCR4BP_2D(t,x,mu,thetaSun0);
        case ['3D' 'CR3BP' 'true']
            rhs = @(t,x) rhs_CR3BP_3D_STM(t,x,mu);
        case ['3D' 'CR3BP' 'false']
            rhs = @(t,x) rhs_CR3BP_3D(t,x,mu);
        case ['3D' 'BCR4BP' 'true']
            rhs = @(t,x) rhs_BCR4BP_3D_STM(t,x,mu,thetaSun0); % STILL TODO
        case ['3D' 'BCR4BP' 'false']
            rhs = @(t,x) rhs_BCR4BP_3D(t,x,mu,thetaSun0);     % STILL TODO
    end
        
    if abs(tf - t0) < 1e-12
        % No integration: the state remains the same
        tt = t0;
        xx = x0(1:n_states)';  
        xf = xx';
        if stm_computation == "true"
            STMf = Phi0;
            STM_hist = Phi0;
        elseif stm_computation == "false"
            STMf = [];
            STM_hist = [];
        end
        return;
    end
    % Perform integration
    options = odeset('reltol', relTol, 'abstol', absTol,'Events',@(t,xx) event_xaxiscrossing(t,xx,evtFlag,evtDir));
    [tt, xx] = ode113(rhs, [t0 tf], x0, options);
    tf = tt(end);

    % Extract state vector and State Transition Matrix
    xf = xx(end,1:n_states)';
    
    if strcmpi(stm_computation,'true')
        STM_hist = rearrange_STM(xx);
        STMf = STM_hist(:,:,end);
    else
        STMf = []; % If STM computation is false, set STMf to empty
        STM_hist = []; % Initialize STM_hist as empty
    end
end

function STM_new = rearrange_STM(xx)
%REARRANGE_STM Recover column-major 4x4 or 6x6 STM slices from augmented solver rows.

    if size(xx,2) == 20
        STM_new = reshape(xx(:,5:end)', 4, 4, []);
    elseif size(xx,2) == 42
        STM_new = reshape(xx(:,7:end)', 6, 6, []);
    else
        error('error: the size of the state vector is neither 20 nor 42');
    end
end
