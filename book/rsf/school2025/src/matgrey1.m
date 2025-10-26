function matgrey1(rsf_file, outfig, figtitle)
    %MATGREY1  Display 2D RSF data in grayscale and save as PDF.
    %   MATGREY1(RSF_FILE, OUTFIG, FIGTITLE) reads RSF data from RSF_FILE,
    %   displays it in grayscale, and saves the figure as a PDF file named OUTFIG.
    %   FIGTITLE is an optional title for the figure.
    %
    % Require 2 or 3 input arguments
    narginchk(2,3);
    if isempty(outfig)
        error('The output filename "outfig" cannot be empty.');
    end
    if nargin < 3
        figtitle = '';
    end
    
    %--- 1. Read RSF header and data
    [data, header] = read_rsf(rsf_file);
    
    taxis = header.o1 + (0:header.n1-1) * header.d1;
    xaxis = header.o2 + (0:header.n2-1) * header.d2;
    
    %--- 2. Display the data in a grayscale figure
    hFig = figure('Color', 'w');
    imagesc(xaxis, taxis, data);
    colormap gray; axis tight;
    xlabel(sprintf('%s (%s)', header.label2, header.unit2));
    ylabel(sprintf('%s (%s)', header.label1, header.unit1)); 
    if ~isempty(figtitle)
        title(figtitle, 'Interpreter', 'none');
    end
    
    %--- 3. Save the figure as a PDF
    if exist('OCTAVE_VERSION','builtin') ~= 0
        set(hFig, 'Units', 'Inches');
        pos = get(hFig, 'Position');  % pos = [left bottom width height]
        set(hFig, ...
            'PaperUnits',   'Inches', ...
            'PaperSize',    [pos(3), pos(4)], ...
            'PaperPosition',[0, 0, pos(3), pos(4)]);
        print(hFig, outfig, '-dpdf', '-painters');
    else
        exportgraphics(hFig, outfig, 'ContentType', 'vector');
    end
    fprintf('Saved figure to %s\n', outfig);
end