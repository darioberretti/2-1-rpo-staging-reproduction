function [xf,STMf,STM_hist, xx, tt, tf]  = propagate(x0,t0,tf,mu,dimensions,framework,stm_computation,varargin)
%PROPAGATE Integrate a state and optional state transition matrix with ode113.
%
% PROTOTYPE:
%   [xf, STMf, STM_hist, xx, tt, tf] = propagate(x0, t0, tf, mu, ...
%       dimensions, framework, stm_computation, evtFlag, evtDir, thetaSun0, ...
%       absTol, relTol)
%   out = propagate(opts)
%
% DESCRIPTION:
%   Integrate planar or spatial CR3BP/BCR4BP motion in the rotating frame.
%   Run01 uses normalized Earth-Moon CR3BP states: positions [DU], velocities
%   [VU], and times [TU]. Automatic STM integration starts from identity and
%   maps initial-state perturbations at t0 to the propagated state.
%
%   The positional interface uses the shared y = 0 crossing event. Struct
%   mode accepts custom events, output times, and an optional custom RHS;
%   it returns exactly one struct containing the six output fields below.
%
% INPUT:
%   x0               [4 or 6 elements] State [x, y, vx, vy] or
%                      [x, y, z, vx, vy, vz]. Converted to a column.
%   t0, tf           Initial and requested final times [TU].
%   mu               Secondary-to-total mass ratio [-].
%   dimensions       Character vector: '2D' or '3D'.
%   framework        Character vector: 'CR3BP' or 'BCR4BP'.
%   stm_computation  Character vector: 'true' or 'false'.
%   evtFlag          Optional logical terminal-event flag; default false.
%   evtDir           Optional crossing direction -1, 0, or 1; default 0.
%   thetaSun0        Optional initial Sun phase [rad]; default 0.
%   absTol, relTol   Optional ode113 tolerances; both default to 1e-12.
%                      Empty positional options select their defaults.
%   opts             Alternative scalar struct with the seven required
%                      fields x0 through stm_computation above. Optional
%                      fields: thetaSun0, absTol, relTol, events, rhs, tspan.
%                      events defaults to [] (no event in struct mode);
%                      tspan defaults to [t0, tf]. rhs is a function handle
%                      @(t, x) with any model parameters already captured.
%                      A custom rhs requires stm_computation = 'false';
%                      supply an augmented state/RHS to integrate a custom STM.
%
% OUTPUT:
%   xf        [n x 1] Final physical state, or full state for a custom RHS.
%   STMf      [n x n] Final STM, or [] when STM integration is disabled.
%   STM_hist  [n x n x Nt] STM at each returned time, or [].
%   xx        [Nt x n] State history, including n^2 appended STM columns
%               when automatic STM integration is enabled.
%   tt        [Nt x 1] Returned integration times [TU].
%   tf        Actual final integration time [TU], possibly event-truncated.
%
%   n is 4 or 6 for built-in dynamics, and Nt is the returned sample count.
%   For abs(tf - t0) < 1e-12, no integration is performed: the state history
%   contains one physical-state row and the STM is identity when requested.

    % Read the struct interface and its optional events, RHS, and output times.
    custom_rhs = [];
    if nargin == 1 && isstruct(x0)
        if nargout > 1
            error('propagate:StructModeSingleOutput', 'Struct-mode propagate() returns a single output struct.');
        end

        opts = x0;
        required = {'x0','t0','tf','mu','dimensions','framework','stm_computation'};
        missing = required(~isfield(opts, required));
        if ~isempty(missing)
            error('propagate:MissingFields', 'Missing required fields: %s', strjoin(missing, ', '));
        end

        x0 = opts.x0;
        t0 = opts.t0;
        tf = opts.tf;
        mu = opts.mu;
        dimensions = opts.dimensions;
        framework = opts.framework;
        stm_computation = opts.stm_computation;

        thetaSun0 = get_opt(opts, 'thetaSun0', 0);
        absTol = get_opt(opts, 'absTol', 1e-12);
        relTol = get_opt(opts, 'relTol', 1e-12);
        eventsFcn = get_opt(opts, 'events', []);
        custom_rhs = get_opt(opts, 'rhs', []);
        tspan_custom = get_opt(opts, 'tspan', []);
        struct_mode = true;
    else
        struct_mode = false;
    end

    % Read positional options, retaining their established defaults.
    if ~struct_mode
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
    end

    % Normalize the state shape and validate custom-RHS compatibility.
    n_states = numel(x0); % number of states
    x0 = x0(:); % verticalize x0

    if ~isempty(custom_rhs) && ~isa(custom_rhs, 'function_handle')
        error('propagate:InvalidCustomRhs', 'opts.rhs must be a function handle.');
    end
    use_custom_rhs = struct_mode && ~isempty(custom_rhs);
    if use_custom_rhs && strcmpi(stm_computation, 'true')
        error('propagate:CustomRhsNoAutoSTM', ...
            ['With custom rhs, set stm_computation=''false'' and pass an ', ...
             'augmented state/rhs if STM integration is needed.']);
    end

    % Append the column-major identity STM when automatic STM integration is enabled.
    if strcmpi(stm_computation,'true')
        Phi0 = eye(n_states);
        x0 = [x0; Phi0(:)];
    end

    % Select the built-in model or use the supplied RHS directly.
    if use_custom_rhs
        rhs = custom_rhs;
    else
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
                rhs = @(t,x) rhs_BCR4BP_3D_STM(t,x,mu,thetaSun0);
            case ['3D' 'BCR4BP' 'false']
                rhs = @(t,x) rhs_BCR4BP_3D(t,x,mu,thetaSun0);
            otherwise
                error('propagate:UnsupportedModel', ...
                    'Unsupported model selection dimensions=%s, framework=%s, stm=%s.', ...
                    dimensions, framework, stm_computation);
        end
    end
        
    % Return the unchanged physical state for a negligible propagation interval.
    if abs(tf - t0) < 1e-12
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
        if struct_mode
            xf = pack_struct_output(xf, STMf, STM_hist, xx, tt, tf);
        end
        return;
    end

    % Configure events/output times, then integrate with the selected tolerances.
    if struct_mode
        if isempty(eventsFcn)
            options = odeset('reltol', relTol, 'abstol', absTol);
        else
            options = odeset('reltol', relTol, 'abstol', absTol, 'Events', eventsFcn);
        end
        if isempty(tspan_custom)
            tspan = [t0 tf];
        else
            tspan = tspan_custom(:).';
            if numel(tspan) < 2
                error('propagate:InvalidTspan', 'opts.tspan must contain at least two time values.');
            end
        end
    else
        options = odeset('reltol', relTol, 'abstol', absTol, ...
            'Events', @(t,xx) event_xaxiscrossing(t,xx,evtFlag,evtDir));
        tspan = [t0 tf];
    end

    [tt, xx] = ode113(rhs, tspan, x0, options);
    tf = tt(end);

    % Recover the final physical state and optional STM history.
    xf = xx(end,1:n_states)';
    
    if strcmpi(stm_computation,'true')
        STM_hist = rearrange_STM(xx);
        STMf = STM_hist(:,:,end);
    else
        STMf = []; % If STM computation is false, set STMf to empty
        STM_hist = []; % Initialize STM_hist as empty
    end

    % Package the single-output struct interface after positional outputs are ready.
    if struct_mode
        xf = pack_struct_output(xf, STMf, STM_hist, xx, tt, tf);
    end
end

function STM_new = rearrange_STM(xx)
%REARRANGE_STM Recover STM slices from augmented ode113 output rows.
%
% xx has 4 + 16 or 6 + 36 columns; STM_new is [n x n x Nt] in
% column-major order, mapping perturbations from the initial epoch.

    if size(xx,2) == 20
        STM_new = reshape(xx(:,5:end)', 4, 4, []);
    elseif size(xx,2) == 42
        STM_new = reshape(xx(:,7:end)', 6, 6, []);
    else
        error('error: the size of the state vector is neither 20 nor 42');
    end
end

function val = get_opt(opts, name, default)
%GET_OPT Use a nonempty struct option or the supplied default.

    if isfield(opts, name) && ~isempty(opts.(name))
        val = opts.(name);
    else
        val = default;
    end
end

function out = pack_struct_output(xf, STMf, STM_hist, xx, tt, tf)
%PACK_STRUCT_OUTPUT Package the six positional outputs for struct-mode callers.

    out = struct('xf', xf, ...
                 'STMf', STMf, ...
                 'STM_hist', STM_hist, ...
                 'xx', xx, ...
                 'tt', tt, ...
                 'tf', tf);
end
