function result = compute_continuation(cfg, map, inputs)
%COMPUTE_CONTINUATION Solve the ascending fixed-altitude warm-start chain.
%
% PROTOTYPE:
%   result = earth2rpo.compute_continuation(cfg, map, inputs)
%
% DESCRIPTION:
%   Build one historical patch seed and solve the configured altitudes in
%   order, passing each final decision vector into the next solve. Retain
%   both stage flags and the node/NLP setup for review. No convergence-based
%   retry or reference comparison is performed inside this campaign.
%
% INPUT:
%   cfg     cfg.earth_to_rpo; continuation.altitudes_km gives the ordered grid
%             and validation supplies the shared solver profile.
%   map     Selected 10000 km Earth-access candidate table.
%   inputs  Tracked planar BCR4BP/CR3BP RPO states in DU/VU.
%
% OUTPUT:
%   result  Struct with data.runs (per-altitude setup/decision vectors/stage
%             flags) and results_table. Table columns use km, m/s, days, and
%             seconds as their suffixes indicate; metrics counts altitudes.

% Require the solver and reconstruct one seed for the ordered altitude campaign.
if ~license('test','Optimization_Toolbox') || isempty(which('fmincon'))
    error('earth2rpo:OptimizationToolbox','Fresh continuation requires Optimization Toolbox.');
end
seed = earth2rpo.build_seed(map,cfg,inputs);
previous = [];
altitudes = cfg.continuation.altitudes_km;
rows = zeros(numel(altitudes),7);

% Solve each altitude and carry its final vector into the next warm start.
for k = 1:numel(altitudes)
    options = cfg.validation;
    options.departure_altitude_km = altitudes(k);
    fprintf('Continuation altitude %.0f km (%d/%d)\n',altitudes(k),k,numel(altitudes));
    solve = solve_vleo_rpo_optimization_2d(options,previous,seed,inputs.x0_rpo_3bp);
    previous = solve.y_opt;
    runs(k) = struct('altitude_km',altitudes(k),'cfg',options, ...
        'constants',solve.setup.constants,'nlp',solve.setup.nlp,'y0',solve.setup.y0, ...
        'stage1_y',solve.stage1.y,'y_opt',solve.y_opt,'exitflag_stage1',solve.stage1.exitflag, ...
        'exitflag_stage2',solve.exitflag_opt); %#ok<AGROW>
    solution = solve.solution;
    rows(k,:) = [altitudes(k),solution.dv_total_mps, ...
        solution.departure_summary.impulse_magnitude_mps,solution.arrival_summary.impulse_magnitude_mps, ...
        solution.tof_days,solution.min_moon_altitude_km,solve.runtime_s];
end

% Store per-altitude solver evidence and a physical-unit comparison table.
result.data.runs = runs;
result.data.results_table = array2table(rows,'VariableNames', ...
    {'altitude_km','total_dv_mps','departure_dv_mps','injection_dv_mps', ...
    'tof_days','min_moon_altitude_km','runtime_s'});
result.metrics.altitude_count = numel(altitudes);
end
