function save_figure_pdf_png(fig, file_base, dpi)
% SAVE_FIGURE_PDF_PNG  Save a figure as PDF and PNG.

if nargin < 3
    dpi = 300;
end
if exist('exportgraphics', 'file') || exist('exportgraphics', 'builtin')
    exportgraphics(fig, [file_base '.pdf'], 'ContentType', 'vector');
    exportgraphics(fig, [file_base '.png'], 'Resolution', dpi);
else
    set(fig, 'Units', 'inches');
    pos = get(fig, 'Position');
    set(fig, 'PaperUnits', 'inches', 'PaperSize', pos(3:4), 'PaperPosition', [0 0 pos(3:4)]);
    print(fig, [file_base '.pdf'], '-dpdf');
    print(fig, [file_base '.png'], '-dpng', sprintf('-r%d', dpi));
end
fprintf('Saved %s.pdf and %s.png\n', file_base, file_base);
end
