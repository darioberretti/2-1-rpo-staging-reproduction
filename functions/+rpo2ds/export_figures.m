function export_figures(figures, output_directory, formats)
%EXPORT_FIGURES Export named runner figures as PNG and optional vector PDF.
%
% PROTOTYPE:
%   rpo2ds.export_figures(figures, output_directory, formats)
%
% DESCRIPTION:
%   Resolve scientific figure numbers with rpo.figure_numbers. Use each
%   figure's configured paper dimensions and close invisible figures after
%   exporting; independent sections retain their global figure numbers.
%
% INPUT:
%   figures           Struct of registered figure names and handles.
%   output_directory  Destination ending in a study folder registered in
%                       rpo.figure_numbers, including 03_earth_to_rpo.
%   formats           png, pdf, both, or a png/pdf string array.
%
% OUTPUT:
%   None. Writes numbered 300 dpi PNG and/or vector PDF files.

% Prepare the destination, requested formats, and fixed study figure numbers.
if ~isfolder(output_directory)
    mkdir(output_directory);
end
formats = lower(string(formats));
figure_numbers = rpo.figure_numbers(output_directory);
if any(formats == "both")
    formats = ["png", "pdf"];
end

% Export at the configured paper size; release invisible batch figures afterwards.
for name_cell = fieldnames(figures).'
    name = name_cell{1};
    handle = figures.(name);
    name = sprintf('%02d_%s', figure_numbers(name), name);
    % Export dimensions are set at creation, independently of the tab size.
    if any(formats == "png")
        % MATLAB's PNG backend can fail while reopening an existing output.
        temporary_png = [tempname(output_directory), '.png'];
        print(handle, temporary_png, '-dpng', '-r300');
        movefile(temporary_png, fullfile(output_directory, name + ".png"), 'f');
    end
    if any(formats == "pdf")
        print(handle, fullfile(output_directory, name + ".pdf"), '-dpdf', '-vector');
    end
    if strcmp(handle.Visible, 'off')
        close(handle);
    end
end
end
