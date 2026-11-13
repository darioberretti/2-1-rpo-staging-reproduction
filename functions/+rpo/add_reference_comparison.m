function result = add_reference_comparison(study, result, reference_filename)
%ADD_REFERENCE_COMPARISON Store and print selected local/reference differences.
%
% PROTOTYPE:
%   result = rpo.add_reference_comparison(study, result, reference_filename)
%
% DESCRIPTION:
%   Compare selected summary metrics against a saved reference result.
%   Differences are descriptive; this function does not apply acceptance
%   tolerances or certify numerical equivalence.
%
% INPUT:
%   study               a1, stability, poincare, or vertical.
%   result              Result struct from the corresponding compute function.
%   reference_filename  MAT file holding the same study result schema.
%
% OUTPUT:
%   result              Input struct with a comparison table containing
%                         quantity, local_value, reference_value, and difference.
%                         Table labels specify the units. Also prints the table.

% Extract matching summary quantities in the same order and units.
reference = rpo.load_result(reference_filename);
switch string(study)
    case "a1"
        labels = ["RPO x0 [DU]"; "RPO vy0 [VU]"; "Period [TU]"; ...
            "Jacobi constant"; "Minimum Moon distance [km]"; ...
            "Closure residual"];
        local_values = local_a1_values(result);
        reference_values = local_a1_values(reference);
    case "stability"
        labels = ["Planar spectral radius"; "Spatial spectral radius"];
        local_values = local_stability_values(result);
        reference_values = local_stability_values(reference);
    case "poincare"
        labels = ["Valid chains"; "Retained crossings"];
        local_values = [result.metrics.valid_chain_count; ...
            numel(result.data.crossing_x)];
        reference_values = [reference.metrics.valid_chain_count; ...
            numel(reference.data.crossing_x)];
    case "vertical"
        labels = ["Maximum |z| [DU]"; "Minimum phase-peak |z| [DU]"; ...
            "Maximum phase-peak |z| [DU]"];
        local_values = local_vertical_values(result);
        reference_values = local_vertical_values(reference);
    otherwise
        error('rpo:ComparisonStudy', 'Unknown comparison study: %s.', study);
end

% Store and print signed differences without applying an acceptance threshold.
comparison = table(labels, local_values, reference_values, ...
    local_values - reference_values, 'VariableNames', ...
    {'quantity', 'local_value', 'reference_value', 'difference'});
result.comparison = comparison;
fprintf('\nLocal/reference comparison\n');
disp(comparison);
end

function values = local_a1_values(result)
%LOCAL_A1_VALUES Collect the six labeled A1/reference-orbit comparison values.

metrics = result.metrics.rpo;
values = [metrics.state(1); metrics.state(4); metrics.period_tu; ...
    metrics.jacobi_constant; metrics.minimum_moon_distance_km; ...
    metrics.closure_residual];
end

function values = local_stability_values(result)
%LOCAL_STABILITY_VALUES Collect planar and spatial reference-orbit spectral radii.

metrics = result.metrics.reference_rpo;
values = [metrics.planar_spectral_radius; metrics.spatial_spectral_radius];
end

function values = local_vertical_values(result)
%LOCAL_VERTICAL_VALUES Collect impulse and phase extrema of absolute z [DU].

impulse_peaks = result.metrics.maximum_absolute_z_by_impulse_du;
phase_peaks = result.metrics.maximum_absolute_z_by_phase_du;
values = [max(impulse_peaks); min(phase_peaks); max(phase_peaks)];
end
