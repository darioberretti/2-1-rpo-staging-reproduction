function stats = statistics(costs,completed)
%STATISTICS Summarize finite completed-trial costs and their empirical distribution.
%
% PROTOTYPE:
%   stats = stationkeeping.statistics(costs, completed)
%
% DESCRIPTION:
%   Filter by completion and finiteness before computing mean,
%   prctile(...,99), and the sorted empirical CDF. Retain all duplicate
%   cost samples. With no valid trials, mean/P99 are NaN and CDF arrays
%   are empty; no distribution fitting or time rescaling is performed.
%
% INPUT:
%   costs      [K x 1] Realized cost per trial [m/s], annual for the
%                frozen 365-day campaign. May include nonfinite entries.
%   completed  [K x 1 logical] True only for trials reaching the horizon.
%
% OUTPUT:
%   stats      Struct with valid_trials, mean_m_s/p99_m_s [m/s],
%                cdf_cost_m_s [M x 1, m/s] and cdf_probability [M x 1].
%                M is the valid-trial count; probabilities are (1:M)/M.

valid = costs(completed & isfinite(costs));
stats = struct('valid_trials',numel(valid),'mean_m_s',NaN,'p99_m_s',NaN, ...
    'cdf_cost_m_s',zeros(0,1),'cdf_probability',zeros(0,1));
if isempty(valid), return; end
stats.mean_m_s = mean(valid);
stats.p99_m_s = prctile(valid,99);
stats.cdf_cost_m_s = sort(valid(:));
stats.cdf_probability = (1:numel(valid)).'/numel(valid);
end
