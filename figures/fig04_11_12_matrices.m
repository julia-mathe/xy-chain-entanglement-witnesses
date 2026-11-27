% FIG04_11_12_MATRICES  Figs. 4, 11 and 12 from data/matrix_figures_N10_T0.1.mat (see fig04_11_12_data).

repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));
outdir = fullfile(repo, 'output', 'fig04_11_12');
if ~exist(outdir, 'dir'), mkdir(outdir); end

MD = load(fullfile(repo, 'data', 'matrix_figures_N10_T0.1.mat'));
plot_matrix_figure(MD.cases([1 2]), {'Ising point', 'XX point'}, fullfile(outdir, 'figC_witness_PBC'));
plot_matrix_figure(MD.cases([3 4]), {'Ising point', 'XX point'}, fullfile(outdir, 'figD_witness_OBC'));
plot_matrix_figure(MD.cases(7), {''}, fullfile(outdir, 'figE_subchain'));
plot_matrix_figure(MD.cases([5 6]), {'', ''}, fullfile(outdir, 'figF_noncritical'));
