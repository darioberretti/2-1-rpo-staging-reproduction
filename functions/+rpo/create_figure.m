function handle = create_figure(visibility,size_pixels,varargin)
%CREATE_FIGURE Keep display tabs together, with an independent export size.
%
% PROTOTYPE:
%   handle = rpo.create_figure(visibility, size_pixels, Name, Value, ...)
%
% DESCRIPTION:
%   Create a docked light-theme figure and set its paper dimensions without
%   setting the display Position. Convert pixel dimensions using the root
%   ScreenPixelsPerInch setting.
%
% INPUT:
%   visibility   Figure visibility, on or off.
%   size_pixels  [1 x 2] Export width/height [pixels], or [] to use the root
%                  DefaultFigurePosition width/height.
%   Name, Value  Additional property pairs forwarded to figure.
%
% OUTPUT:
%   handle       Figure handle with paper dimensions in inches.

% Resolve export dimensions independently of the docked figure tab.
if isempty(size_pixels)
    default_position = get(groot,'DefaultFigurePosition');
    size_pixels = default_position(3:4);
end
size_inches = size_pixels / get(groot,'ScreenPixelsPerInch');
% Setting figure Position would undock it from MATLAB's figure container.
handle = figure('Visible',visibility,'Theme','light','Color','w', ...
    'WindowStyle','docked','PaperUnits','inches', ...
    'PaperPosition',[0,0,size_inches],'PaperSize',size_inches, ...
    'PaperPositionMode','manual',varargin{:});
end
