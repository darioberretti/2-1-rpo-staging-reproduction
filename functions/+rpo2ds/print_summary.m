function print_summary(result, section_name)
%PRINT_SUMMARY Print the compact Run 02 section metrics.
%
% PROTOTYPE:
%   rpo2ds.print_summary(result, section_name)
%
% DESCRIPTION:
%   Print stored metrics without recomputing or validating them. Sections
%   without a metrics field return without printing a summary.
%
% INPUT:
%   result        Compact validation, flexibility, or reachability result.
%   section_name  venus, mars, flexibility, or reachability.
%
% OUTPUT:
%   None. Reports delta-v [m/s], flight time [days], speeds [km/s], and/or
%   failed case indices as applicable.

switch lower(string(section_name))
    case {"venus", "mars"}
        if ~isfield(result, 'metrics'), return, end
        m = result.metrics;
        fprintf('%s validation: total %.3f m/s, RPO %.3f m/s, Earth %.3f m/s, TOF %.3f d.\n', ...
            char(upper(string(section_name))), m.total_delta_v_m_s, m.rpo_departure_delta_v_m_s, ...
            m.earth_perigee_delta_v_m_s, m.time_of_flight_days);
    case "flexibility"
        if ~isfield(result, 'metrics'), return, end
        for name = ["venus", "mars"]
            m = result.metrics.(name);
            fprintf('%s flexibility: %d converged, dv [%.3f, %.3f] m/s, TOF [%.3f, %.3f] d, failed indices %s.\n', ...
                upper(name), m.converged_count, m.delta_v_range_m_s, ...
                m.time_of_flight_range_days, local_indices(m.failed_indices));
        end
    case "reachability"
        if ~isfield(result, 'metrics'), return, end
        fprintf('Reachability: |v_inf| [%.6f, %.6f] km/s; failed indices %s.\n', ...
            result.metrics.v_inf_range_km_s, local_indices(result.metrics.failed_indices));
end
end

function value = local_indices(indices)
%LOCAL_INDICES Format failed case indices as a row-vector string, or [].

if isempty(indices)
    value = '[]';
else
    value = mat2str(indices(:).');
end
end
