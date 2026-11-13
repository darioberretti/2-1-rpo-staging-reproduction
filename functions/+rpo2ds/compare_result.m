function comparison = compare_result(candidate, reference_filename, section_name, destination)
%COMPARE_RESULT Compare compact candidate/reference numerical payloads.
%
% PROTOTYPE:
%   comparison = rpo2ds.compare_result(candidate, reference_filename, ...
%       section_name, destination)
%
% DESCRIPTION:
%   Recursively compare shared numeric/logical fields, nested structs, and
%   matching numeric cells. Report shapes and maximum absolute/relative
%   differences without applying acceptance tolerances. Nonshared fields
%   and text are not compared; a missing reference skips the comparison.
%
% INPUT:
%   candidate           Candidate result struct.
%   reference_filename  MAT file holding the reference result.
%   section_name        Label for console output.
%   destination         Optional comparison MAT path; [] disables saving.
%
% OUTPUT:
%   comparison          Struct with reference_available, per-field paths/shape
%                         flags/differences, field count, and aggregate maxima.
%                         Absolute differences retain each field's own units;
%                         relative differences are scaled by max(abs(ref), eps).

comparison = struct('section', string(section_name), 'reference_available', false, ...
    'fields', struct('path', {}, 'same_shape', {}, 'max_abs', {}, 'max_rel', {}), ...
    'compared_field_count', 0, 'maximum_absolute_difference', NaN, ...
    'maximum_relative_difference', NaN);
if ~isfile(reference_filename)
    fprintf('No accepted reference result is available for %s; comparison skipped.\n', section_name);
    return
end
reference = rpo2ds.load_result(reference_filename);
comparison.reference_available = true;
comparison.fields = local_compare_struct(candidate, reference, 'result');
comparison.compared_field_count = numel(comparison.fields);
if ~isempty(comparison.fields)
    comparison.maximum_absolute_difference = max([comparison.fields.max_abs], [], 'omitnan');
    comparison.maximum_relative_difference = max([comparison.fields.max_rel], [], 'omitnan');
end
fprintf('%s comparison: %d numeric fields, max abs %.3e, max rel %.3e.\n', ...
    section_name, comparison.compared_field_count, ...
    comparison.maximum_absolute_difference, comparison.maximum_relative_difference);
if nargin >= 4 && ~isempty(destination)
    parent = fileparts(destination);
    if ~isfolder(parent), mkdir(parent); end
    save(destination, 'comparison');
end
end

function fields = local_compare_struct(candidate, reference, prefix)
%LOCAL_COMPARE_STRUCT Walk shared fields and collect numeric difference summaries.
%
% PROTOTYPE:
%   fields = local_compare_struct(candidate, reference, prefix)
%
% DESCRIPTION:
%   Compare matching struct fields recursively, including equal-size struct
%   arrays and numeric cells. Record numeric shape mismatches; skip text,
%   nonshared fields, and unsupported cell contents.
%
% INPUT:
%   candidate, reference  Scalar structs at the current recursion level.
%   prefix                Text path identifying the containing field/element.
%
% OUTPUT:
%   fields                Struct array with path, same_shape, max_abs, max_rel.
%                         Differences use native field units and the
%                         reference-value relative scaling.

fields = struct('path', {}, 'same_shape', {}, 'max_abs', {}, 'max_rel', {});
if ~isstruct(candidate) || ~isstruct(reference)
    return
end
names = intersect(fieldnames(candidate), fieldnames(reference), 'stable');
for index = 1:numel(names)
    name = names{index};
    left = candidate.(name);
    right = reference.(name);
    path = string(prefix) + "." + string(name);
    if (isnumeric(left) || islogical(left)) && ...
            (isnumeric(right) || islogical(right))
        entry = struct('path', path, 'same_shape', isequal(size(left), size(right)), ...
            'max_abs', NaN, 'max_rel', NaN);
        if entry.same_shape && ~isempty(left)
            difference = abs(double(left) - double(right));
            scale = max(abs(double(right(:))), eps);
            entry.max_abs = max(difference(:));
            entry.max_rel = max((difference(:) ./ scale(:)));
        end
        fields(end + 1) = entry; %#ok<AGROW>
    elseif isstruct(left) && isstruct(right)
        if ~isequal(size(left), size(right))
            fields(end + 1) = struct('path', path, 'same_shape', false, ...
                'max_abs', NaN, 'max_rel', NaN); %#ok<AGROW>
        else
            for element_index = 1:numel(left)
                element_path = path;
                if ~isscalar(left)
                    element_path = path + "(" + element_index + ")";
                end
                nested = local_compare_struct(left(element_index), ...
                    right(element_index), element_path);
                fields = [fields(:); nested(:)];
            end
        end
    elseif iscell(left) && iscell(right) && numel(left) == numel(right)
        for cell_index = 1:numel(left)
            if (isnumeric(left{cell_index}) || islogical(left{cell_index})) && ...
                    (isnumeric(right{cell_index}) || islogical(right{cell_index})) && ...
                    isequal(size(left{cell_index}), size(right{cell_index}))
                difference = abs(double(left{cell_index}) - double(right{cell_index}));
                scale = max(abs(double(right{cell_index}(:))), eps);
                fields(end + 1) = struct('path', path + "{" + cell_index + "}", ...
                    'same_shape', true, 'max_abs', max(difference(:)), ...
                    'max_rel', max(difference(:) ./ scale(:))); %#ok<AGROW>
            end
        end
    end
end
end
