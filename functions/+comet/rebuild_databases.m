function rebuild_databases(settings, directory)
%REBUILD_DATABASES Compute and save the heliocentric cumulative-cost database.
%
% PROTOTYPE:
%   comet.rebuild_databases(settings, directory)
%
% DESCRIPTION:
%   Build all three departure grids, then save the completed database.
%   Existing output is replaced only after the new MAT file saves.
%   comet.prepare_run enforces the user's replacement policy beforehand.
%
% INPUT:
%   settings   cfg.comet_case_study from reproduction_config().
%              heliocentric and constants define the frozen grids/model.
%   directory  Destination folder for comet.database_file.
%
% OUTPUT:
%   None. Writes one v7.3 MAT file containing the database struct whose
%   fields, dimensions and units are documented in comet.build_database.

if ~isfolder(directory), mkdir(directory); end
database = comet.build_database(settings.heliocentric, settings.constants);
filename = fullfile(directory, comet.database_file());
temporary = [tempname(directory), '.mat'];
save(temporary, 'database', '-v7.3');
movefile(temporary, filename, 'f');
fprintf('Heliocentric database saved in %s\n', filename);
end
