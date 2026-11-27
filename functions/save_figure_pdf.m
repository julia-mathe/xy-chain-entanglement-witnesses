function save_figure_pdf(figure_handle, file_name, width_inches, height_inches)
% SAVE_FIGURE_PDF  Save a figure as a PDF of given size (inches).

set(figure_handle, 'PaperUnits', 'inches', ...
                   'PaperSize', [width_inches, height_inches], ...
                   'PaperPosition', [0, 0, width_inches, height_inches]);

if exist('OCTAVE_VERSION', 'builtin')
    print(figure_handle, file_name, '-dpdf', '-nosvgconvert');
else
    print(figure_handle, file_name, '-dpdf');
end
fprintf('Saved %s\n', file_name);
end
