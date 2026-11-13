function samples_years = comet_sample_warning_times(law, n_samples)
%COMET_SAMPLE_WARNING_TIMES Sample a piecewise-uniform warning-time law without reordering RNG draws.
%
% PROTOTYPE:
%   samples_years = comet_sample_warning_times(law, n_samples)
%
% DESCRIPTION:
%   For each sample, draw the probability segment first and then
%   its uniform position. The caller owns RNG initialization; retain
%   this serial two-draw sequence for reproducible shared catalogues.
%
% INPUT:
%   law        Struct with segment_edges_years [1 x K+1, years] and
%                segment_masses [1 x K, probabilities summing to one].
%   n_samples  Number of requested warning samples.
%
% OUTPUT:
%   samples_years  [n_samples x 1] Warning times [years].

    edges = law.segment_edges_years;                                   % segment edges [years]
    masses = law.segment_masses;                                       % segment masses
    cdf = cumsum(masses);                                              % cumulative distribution of the segment masses
    samples_years = zeros(n_samples, 1);                               % sampled warning times [years]

    % Draw the segment index and then sample uniformly inside the segment.
    for i_sample = 1:n_samples
        r = rand();                                                    % segment-level random number
        idx = find(r <= cdf, 1, 'first');                              % selected segment index
        left_edge = edges(idx);                                        % selected segment left edge [years]
        right_edge = edges(idx + 1);                                   % selected segment right edge [years]
        samples_years(i_sample) = left_edge + (right_edge - left_edge) * rand();
    end
end
