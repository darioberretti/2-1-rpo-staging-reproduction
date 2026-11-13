function available = parallel_available(requested)
%PARALLEL_AVAILABLE True when parallel execution is requested and installed.
%
% requested is a scalar logical. available also requires the Parallel
% Computing Toolbox license and installation; this does not start a pool.

available = requested && license('test', 'Distrib_Computing_Toolbox') ...
    && ~isempty(ver('parallel'));
end
