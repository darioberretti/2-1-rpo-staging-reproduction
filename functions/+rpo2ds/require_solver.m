function require_solver(backend)
%REQUIRE_SOLVER Fail before an expensive section when its solver is absent.
%
% PROTOTYPE:
%   rpo2ds.require_solver(backend)
%
% DESCRIPTION:
%   Check the configured backend. For CasADi/IPOPT, try the repository-local
%   third_party/casadi installation if casadi.Opti is not already on the path.
%   This presence check does not execute an optimization or test convergence.
%
% INPUT:
%   backend  casadi_ipopt or fmincon.
%
% OUTPUT:
%   None. May add a local dependency path; errors when the backend is absent.

backend = lower(string(backend));
switch backend
    case "casadi_ipopt"
        if isempty(which('casadi.Opti'))
            repository_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            local_casadi = fullfile(repository_root, 'third_party', 'casadi');
            if isfolder(local_casadi)
                addpath(genpath(local_casadi), '-begin');
            end
        end
        if isempty(which('casadi.Opti'))
            error('rpo2ds:MissingCasadi', [ ...
                'CasADi/IPOPT is required for fresh transfer computation, ', ...
                'but casadi.Opti is not available. Install CasADi and add its ', ...
                'MATLAB interface to the path or place it under third_party/casadi ', ...
                'before rerunning. Saved-result ', ...
                'replay does not require CasADi.']);
        end
    case "fmincon"
        if isempty(which('fmincon'))
            error('rpo2ds:MissingOptimizationToolbox', [ ...
                'Optimization Toolbox (fmincon) is required for this fresh ', ...
                'transfer computation. Install or enable it before rerunning. ', ...
                'Saved-result replay does not require it.']);
        end
    otherwise
        error('rpo2ds:UnknownSolverBackend', 'Unsupported Run 02 solver backend: %s.', backend);
end
end
