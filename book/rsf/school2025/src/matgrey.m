function matgrey(rsf_file, outfig, figtitle)
    % Require 2 or 3 input arguments
    narginchk(2,3);
    if isempty(outfig)
            error('The output filename "outfig" cannot be empty.');
    end
    if nargin < 3
            figtitle = '';
    end
    %--- 1. Read RSF header and data
    % Add the RSF library to MATLAB path
    addpath(fullfile(getenv('RSFROOT'), 'lib'));
    [shape, delta, origin, label, unit] = rsf_read_header(rsf_file);
    data = zeros(shape);
    rsf_read(data, rsf_file);
    time_axis = origin(1) + (0:shape(1)-1) * delta(1);
    spatial_axis = origin(2) + (0:shape(2)-1) * delta(2);
    %--- 2. Display the data in a grayscale figure
    hFig = figure('Color', 'w');
    imagesc(spatial_axis, time_axis, data);
    colormap gray;        axis tight;
    xlabel(sprintf('%s (%s)', label{2}, unit{2}));
    ylabel(sprintf('%s (%s)', label{1}, unit{1})); 
    if ~isempty(figtitle)
            title(figtitle, 'Interpreter', 'none');
    end
    %--- 3. Save the figure as a PDF
    exportgraphics(hFig, outfig, 'ContentType', 'vector');
    fprintf('Saved figure to %s\n', outfig);
end
