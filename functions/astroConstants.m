function value = astroConstants(identifier)
%ASTROCONSTANTS Return the historical constants retained by the Run 02 backend.
%
% PROTOTYPE:
%   value = astroConstants(identifier)
%
% DESCRIPTION:
%   Provide only the constant identifiers used by the migrated workflow;
%   preserve its historical numerical values without an ephemeris library.
%
% INPUT:
%   identifier  Scalar/vector: 4 Sun GM, 13 Earth GM, 20 Moon GM,
%                23 Earth mean radius, or 30 Moon mean radius.
%
% OUTPUT:
%   value       Row vector in input traversal order. GM is in km^3/s^2
%                and radii in km. Unsupported identifiers raise an error.

identifiers = identifier(:).';
value = zeros(size(identifiers));
for index = 1:numel(identifiers)
    switch identifiers(index)
        case 4
            value(index) = 1.32712440017987e11; % Sun GM [km^3/s^2]
        case 13
            value(index) = 3.98600433e5;         % Earth GM [km^3/s^2]
        case 20
            value(index) = 4902.801;              % Moon GM [km^3/s^2]
        case 23
            value(index) = 6371.01;               % Earth mean radius [km]
        case 30
            value(index) = 1738.0;                % Moon mean radius [km]
        otherwise
            error('rpo2ds:UnsupportedConstant', ...
                'astroConstants identifier %d is not part of the public Run 02 backend.', ...
                identifiers(index));
    end
end
end
