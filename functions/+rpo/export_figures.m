function export_figures(figures, output_directory, formats)
%EXPORT_FIGURES Export named figures as publication-quality PNG and/or PDF.
%
% PROTOTYPE:
%   rpo.export_figures(figures, output_directory, formats)
%
% DESCRIPTION:
%   Use the study-wide numbering from rpo.figure_numbers, independent of
%   which subset of figures is supplied. Export at the figures' paper sizes
%   as vector PDF and/or 300 dpi PNG; close invisible figures after export.
%
% INPUT:
%   figures           Struct mapping registered figure names to handles.
%   output_directory  Destination path ending in a registered study folder.
%   formats           png, pdf, both, or an array containing png and pdf.
%
% OUTPUT:
%   None. Creates numbered files in output_directory.

% Prepare the destination and resolve the complete study-wide figure numbering.
if ~isfolder(output_directory)
    mkdir(output_directory);
end
names = fieldnames(figures);
figure_numbers = rpo.figure_numbers(output_directory);
formats = lower(string(formats));
if any(formats == "both")
    formats = ["png", "pdf"];
end

% Export each supplied figure under its fixed number and requested formats.
for index = 1:numel(names)
    name = names{index};
    figure_handle = figures.(name);
    name = sprintf('%02d_%s', figure_numbers(name), name);
    % Export dimensions are set at creation, independently of the tab size.
    if any(formats == "pdf")
        print(figure_handle, fullfile(output_directory, ...
            string(name) + ".pdf"), '-dpdf', '-vector');
    end
    if any(formats == "png")
        print(figure_handle, fullfile(output_directory, ...
            string(name) + ".png"), '-dpng', '-r300');
    end

    % Release batch figures while leaving visible interactive tabs open.
    if strcmp(figure_handle.Visible, 'off')
        close(figure_handle);
    end
end
end
