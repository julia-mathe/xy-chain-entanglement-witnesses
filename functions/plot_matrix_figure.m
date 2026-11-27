function fig = plot_matrix_figure(records, row_labels, file_base)
% PLOT_MATRIX_FIGURE  Draws Figs. 4, 11 and 12: Gamma, Z* and the distance profiles of Z*.

nrows = numel(records);
FONT = 'Times New Roman';
colors  = [ 42 120 214;
           235 104  52;
            27 175 122] / 255;
markers = {'o', 's', '^'};
dodge   = [-0.17, 0, 0.17];

width = 7.4;  row_h = 2.25;  top = 0.12;
s = 1.52;
x_axes = [0.55, 3.02, 5.84];
height = top + nrows * row_h;
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 width height]);
cmap = rdbu_colormap(256);

for i = 1:nrows
    Gamma = records(i).Gamma;
    Z = records(i).Z;
    N = size(Gamma, 1) / 3;
    perm = [1:3:3*N, 2:3:3*N, 3:3:3*N];
    y0 = height - top - i*row_h + 0.48;

    titles = {'$\Gamma$', sprintf('$Z^\\star,\\ \\omega_{\\rm CMC}=%.3f$', records(i).omega)};
    mats = {Gamma(perm, perm), Z(perm, perm)};
    for p = 1:2
        ax = axes('Parent', fig, 'Units', 'inches', 'Position', [x_axes(p), y0, s, s]);
        M = mats{p};
        v = max(abs(M(:)));
        if v < 1e-12, v = 1; end
        imagesc(ax, M);
        colormap(ax, cmap);
        caxis(ax, [-v, v]);
        axis(ax, 'image');
        hold(ax, 'on');
        for kk = 1:2
            plot(ax, [kk*N + 0.5, kk*N + 0.5], [0.5, 3*N + 0.5], 'k-', 'LineWidth', 0.8);
            plot(ax, [0.5, 3*N + 0.5], [kk*N + 0.5, kk*N + 0.5], 'k-', 'LineWidth', 0.8);
        end
        centers = [N/2, 3*N/2, 5*N/2] + 0.5;
        set(ax, 'XTick', centers, 'YTick', centers, 'XTickLabel', {'$x$', '$y$', '$z$'}, ...
            'YTickLabel', {'$x$', '$y$', '$z$'}, 'TickLabelInterpreter', 'latex', ...
            'TickLength', [0 0], 'FontName', FONT, 'FontSize', 12, 'LineWidth', 0.8, 'Box', 'on');
        title(ax, titles{p}, 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'normal');
        pos = get(ax, 'Position');
        cb = colorbar(ax);
        set(ax, 'Position', pos);
        set(cb, 'Units', 'inches', 'Position', [x_axes(p) + s + 0.06, y0, 0.09, s]);
        ticks = nice_ticks(-v, v);
        if isprop(cb, 'Ticks')
            set(cb, 'Ticks', ticks);
        else
            set(cb, 'YTick', ticks);
        end
        set(cb, 'FontName', FONT, 'FontSize', 9, 'LineWidth', 0.6);
        if p == 1 && ~isempty(row_labels{i})
            ylabel(ax, row_labels{i}, 'FontName', FONT, 'FontWeight', 'bold', 'FontSize', 13);
        end
    end

    ax = axes('Parent', fig, 'Units', 'inches', 'Position', [x_axes(3), y0, s, s]);
    hold(ax, 'on');
    r = 0:N-1;
    lo = 0;  hi = 0;
    for alpha = 1:3
        B = Z(alpha:3:3*N, alpha:3:3*N);
        profile = zeros(1, N);
        for d = r
            profile(d + 1) = mean([diag(B, d); diag(B, -d)]);
        end
        lo = min([lo, profile]);  hi = max([hi, profile]);
        plot(ax, r + dodge(alpha), profile, '-', 'Color', colors(alpha, :), 'Marker', markers{alpha}, ...
             'MarkerSize', 4.5, 'MarkerFaceColor', colors(alpha, :), 'LineWidth', 1.0);
    end
    plot(ax, [-0.5, N - 0.5], [0, 0], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.7);
    uistack_bottom(ax);
    pad = 0.06 * max(hi - lo, 1e-3);
    ylim_ = [lo - pad, hi + pad];
    set(ax, 'XLim', [-0.5, N - 0.5], 'XTick', r, 'YLim', ylim_, 'YTick', nice_ticks(ylim_(1), ylim_(2)), ...
        'FontName', FONT, 'FontSize', 10, 'LineWidth', 0.8, 'Box', 'on', 'TickDir', 'in');
    grid(ax, 'on');
    set(ax, 'GridAlpha', 0.25);
    xlabel(ax, 'distance $r$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel(ax, '$\overline{(Z^\star)^{\alpha\alpha}}(r)$', 'Interpreter', 'latex', 'FontSize', 12);
    if i == 1
        draw_legend(fig, [x_axes(3), y0 + s + 0.03, s, 0.2], colors, markers, FONT);
    end
end

if nargin >= 3 && ~isempty(file_base)
    save_figure_pdf_png(fig, file_base, 400);
end
end

function ticks = nice_ticks(lo, hi)
range = hi - lo;
e = 10^floor(log10(range / 4));
steps = [1, 2, 5, 10, 20, 50] * e;
k = find(range ./ steps <= 4 + 1e-9, 1);
ticks = ticks_with_step(lo, hi, steps(k));
if numel(ticks) < 3 && k > 1
    ticks = ticks_with_step(lo, hi, steps(k - 1));
end
end

function ticks = ticks_with_step(lo, hi, step)
ticks = step * (ceil(lo/step - 1e-9) : floor(hi/step + 1e-9));
ticks(abs(ticks) < 1e-9 * step) = 0;
end

function draw_legend(fig, position, colors, markers, font)
ax = axes('Parent', fig, 'Units', 'inches', 'Position', position, 'Visible', 'off', ...
          'XLim', [0 1], 'YLim', [0 1]);
hold(ax, 'on');
labels = {'$x$', '$y$', '$z$'};
for alpha = 1:3
    x = 0.05 + 0.33 * (alpha - 1);
    plot(ax, [x, x + 0.16], [0.5, 0.5], '-', 'Color', colors(alpha, :), 'LineWidth', 1.0);
    plot(ax, x + 0.08, 0.5, markers{alpha}, 'Color', colors(alpha, :), ...
         'MarkerFaceColor', colors(alpha, :), 'MarkerSize', 4.5);
    text(ax, x + 0.20, 0.5, labels{alpha}, 'Interpreter', 'latex', 'FontName', font, ...
         'FontSize', 12, 'VerticalAlignment', 'middle');
end
end

function uistack_bottom(ax)
children = get(ax, 'Children');
if numel(children) > 1
    set(ax, 'Children', [children(2:end); children(1)]);
end
end
