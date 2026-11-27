% FIG05_STRUCTURE_FACTORS  Fig. 5: structure factors at T = 0.1, infinite chain and N = 10.
% No CVX needed.

%% Setup
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));
if ~exist(fullfile(repo, 'output'), 'dir'), mkdir(fullfile(repo, 'output')); end

T = 0.1;
colors = [ 31  90 166
          209  73  91
           42 157 143] / 255;

panels = {0, 0, 'XX point (\gamma, h/J) = (0, 0)',    [pi, pi, 0]
          1, 1, 'Ising point (\gamma, h/J) = (1, 1)', [pi, 0, pi]};

figure_handle = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 6.3, 2.55]);

for n = 1:2
    gamma = panels{n, 1};
    h = panels{n, 2};
    selected_q = panels{n, 4};

    [mz, Cxx, Cyy, Czz] = spin_correlations_tdl(700, gamma, h, T);
    Czc = Czz - mz^2;
    Czc(1) = 1 - mz^2;
    C = {Cxx, Cyy, Czc};
    q = linspace(0, pi, 601);

    ring = ring_spectrum(10, gamma, h);
    [~, ~, ~, S10, q10] = ring_cmc_at_temperature(ring, T);
    keep = q10 <= pi + 1e-9;

    panel_left = [0.085, 0.545];
    axes('Position', [panel_left(n), 0.19, 0.43, 0.70]);
    hold on;
    box on;
    line_handles = zeros(1, 3);
    for alpha = 1:3
        S = structure_factor_from_correlations(C{alpha}, q);
        S(end) = staggered_structure_factor_euler(C{alpha});
        if alpha == 2 && gamma == 0
            style = '--';
        else
            style = '-';
        end
        line_handles(alpha) = plot(q/pi, S, style, 'Color', colors(alpha, :), 'LineWidth', 1.3);
        plot(q10(keep)/pi, S10(alpha, keep), 'o', 'MarkerSize', 3.5, 'MarkerFaceColor', 'w', ...
             'MarkerEdgeColor', colors(alpha, :), 'LineWidth', 0.9);
        qa = selected_q(alpha);
        if qa < pi
            Sa = interp1(q, S, qa);
        else
            Sa = S(end);
        end
        plot(qa/pi, Sa, 'p', 'MarkerSize', 9, 'MarkerFaceColor', colors(alpha, :), ...
             'MarkerEdgeColor', 'k', 'LineWidth', 0.4);
    end
    set(gca, 'YScale', 'log', 'XLim', [-0.015, 1.015], 'YLim', [0.02, 40], ...
             'XTick', 0:0.25:1, 'FontName', 'Times', 'FontSize', 8);
    xlabel('q/\pi', 'FontSize', 9);
    title(panels{n, 3}, 'FontWeight', 'normal', 'FontSize', 9);
    if n == 2
        set(gca, 'YTickLabel', []);
    end
    if n == 1
        ylabel('structure factor S_\alpha(q)', 'FontSize', 9);
        legend_handle = legend(line_handles, {'S_x(q)', 'S_y(q)', 'S_z(q)'}, ...
                               'Location', 'north', 'Orientation', 'horizontal');
        set(legend_handle, 'Box', 'off', 'FontSize', 8);
    end
end

save_figure_pdf(figure_handle, fullfile(repo, 'output', 'structure_factors_T0.1.pdf'), 6.3, 2.55);
