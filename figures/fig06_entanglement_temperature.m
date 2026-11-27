% FIG06_ENTANGLEMENT_TEMPERATURE  Fig. 6 from the data written by the fig06_data_* scripts.

%% Setup
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));
if ~exist(fullfile(repo, 'output'), 'dir'), mkdir(fullfile(repo, 'output')); end

names = {'entanglement_temperature_tdl.mat', 'entanglement_temperature_N10.mat', ...
         'entanglement_temperature_N4.mat', 'energy_witness_threshold.mat', 'omega_vs_temperature.mat'};
for n = 1:numel(names)
    if ~exist(fullfile(repo, 'data', names{n}), 'file')
        error('Data file data/%s not found: run the fig06_data_* scripts first.', names{n});
    end
end
tdl = load(fullfile(repo, 'data', names{1}));
n10 = load(fullfile(repo, 'data', names{2}));
n4  = load(fullfile(repo, 'data', names{3}));
ew  = load(fullfile(repo, 'data', names{4}));
ot  = load(fullfile(repo, 'data', names{5}));

blue = [31 90 166] / 255;
red  = [209 73 91] / 255;
line_colors = {blue, red};

figure_handle = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 6.3, 2.6]);

%% (a) entanglement temperature
axes_a = axes('Position', [0.08, 0.17, 0.40, 0.78]);
hold on;
box on;
plot([1, 1], [0, 2.8], '-', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 0.6);
color_handles = zeros(1, 2);
for ig = 1:2
    c = line_colors{ig};
    color_handles(ig) = plot(tdl.fields, tdl.T_ent(ig, :), '-', 'Color', c, 'LineWidth', 1.3);
    plot(n10.fields, n10.T_ent(ig, :), '--', 'Color', c, 'LineWidth', 1.0);
    plot(n4.fields, n4.T_ent(ig, :), ':', 'Color', c, 'LineWidth', 1.2);
    plot(ew.fields, ew.T_E(ig, :), '-.', 'Color', c, 'LineWidth', 0.8);
end
set(axes_a, 'XLim', [0, 2.5], 'YLim', [0, 2.8], 'FontName', 'Times', 'FontSize', 8);
xlabel('h/J', 'FontSize', 9);
ylabel('T_{ent}/J', 'FontSize', 9);
text(0.03, 0.93, '(a)', 'Units', 'normalized', 'FontName', 'Times');
legend_colors = legend(color_handles, {'XX line, \gamma = 0', 'Ising line, \gamma = 1'}, 'Location', 'northwest');
set(legend_colors, 'Box', 'off', 'FontSize', 7.5);

axes_styles = axes('Position', get(axes_a, 'Position'), 'Visible', 'off');
hold on;
gray = [0.3, 0.3, 0.3];
style_handles = [plot(nan, nan, '-', 'Color', gray, 'LineWidth', 1.3), ...
                 plot(nan, nan, '--', 'Color', gray, 'LineWidth', 1.0), ...
                 plot(nan, nan, ':', 'Color', gray, 'LineWidth', 1.2), ...
                 plot(nan, nan, '-.', 'Color', gray, 'LineWidth', 0.8)];
legend_styles = legend(style_handles, {'CMC, N = \infty', 'CMC, N = 10', 'CMC, N = 4', 'energy witness'}, ...
                       'Location', 'southeast');
set(legend_styles, 'Box', 'off', 'FontSize', 7);

%% (b) omega_CMC as a function of temperature
axes_b = axes('Position', [0.58, 0.17, 0.40, 0.78]);
hold on;
box on;
fill([0, 2.6, 2.6, 0], [1, 1, 1.25, 1.25], [0.9, 0.9, 0.9], 'EdgeColor', 'none');
plot([0, 2.6], [1, 1], 'k-', 'LineWidth', 0.6);
h1 = plot(ot.temperatures, ot.omega_infinite(1, :), '-', 'Color', blue, 'LineWidth', 1.3);
h2 = plot(ot.temperatures, ot.omega_N10(1, :), '--', 'Color', blue, 'LineWidth', 1.0);
h3 = plot(ot.temperatures, ot.omega_infinite(2, :), '-', 'Color', red, 'LineWidth', 1.3);
h4 = plot(ot.temperatures, ot.omega_N10(2, :), '--', 'Color', red, 'LineWidth', 1.0);
text(0.12, 1.08, 'no detection', 'FontSize', 8, 'Color', [0.35, 0.35, 0.35], 'FontName', 'Times');
set(axes_b, 'XLim', [0, 2.6], 'YLim', [0, 1.25], 'FontName', 'Times', 'FontSize', 8);
xlabel('T/J', 'FontSize', 9);
ylabel('\omega_{CMC}', 'FontSize', 9);
text(0.03, 0.93, '(b)', 'Units', 'normalized', 'FontName', 'Times');
legend_b = legend([h1, h2, h3, h4], {'XX point, N = \infty', 'N = 10', 'Ising point, N = \infty', 'N = 10'}, ...
                  'Location', 'southeast');
set(legend_b, 'Box', 'off', 'FontSize', 7.5);

save_figure_pdf(figure_handle, fullfile(repo, 'output', 'entanglement_temperature.pdf'), 6.3, 2.6);
