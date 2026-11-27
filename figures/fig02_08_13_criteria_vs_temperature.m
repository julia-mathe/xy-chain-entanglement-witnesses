% FIG02_08_13_CRITERIA_VS_TEMPERATURE  Figs. 2, 8 and 13 from data/N4_criteria_vs_T_<BC>.mat
% (see fig02_08_13_data). MATLAB R2020b or newer.

repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));
outdir = fullfile(repo, 'output', 'fig02_08_13');
if ~exist(outdir, 'dir'), mkdir(outdir); end

D = cell(1, 2);
for ib = 1:2
    bc = {'open', 'periodic'};
    D{ib} = load_criteria_data(fullfile(repo, 'data', sprintf('N4_criteria_vs_T_%s.mat', bc{ib})));
end

plot_criteria(D, [0 0; 1 1], {'XX point, $(\gamma,h)=(0,0)$', 'Ising point, $(\gamma,h)=(1,1)$'}, ...
    fullfile(outdir, 'criteria_strengths_vs_T_N4'));
plot_criteria(D, [0.6 0.2; 0.9 0.9; 1.3 1.2], {'$(\gamma,h)=(0.6,\,0.2)$', '$(\gamma,h)=(0.9,\,0.9)$', ...
    '$(\gamma,h)=(1.3,\,1.2)$'}, fullfile(outdir, 'criteria_strengths_vs_T_N4_appendix'));
plot_gme(D, [0 0; 0.6 0.2; 0.9 0.9; 1.3 1.2], fullfile(outdir, 'gme_ppt_mixer_vs_T_N4'));

function D = load_criteria_data(file)
S = load(file, 'settings', 'results');
R = S.results(1);
D = struct('ENTSDP', R.ENTSDP, 'ENTCMC', R.ENTCMC, 'ENTPPT', R.ENTPPT, ...
    'FPPTwitStrength', R.FPPTwitStrength, 'FPPTwitRaw', R.FPPTwitRaw, ...
    'T_list', S.settings.T_list, 'parameterCases', S.settings.parameterCases, ...
    'detectTol', S.settings.detectTol);
end

function ic = case_index(D, gh)
ic = find(all(abs(D.parameterCases - gh) < 1e-9, 2), 1);
if isempty(ic)
    error('No data for (gamma, h) = (%g, %g); run fig02_08_13_data first.', gh(1), gh(2));
end
end

function [y, Tlast] = detected_series(D, key, ic)
y = D.(key)(ic, :);
if strcmp(key, 'FPPTwitStrength')
    y(~isfinite(D.FPPTwitRaw(ic, :))) = NaN;
end
det = isfinite(y) & (y > D.detectTol);
y(~det) = NaN;
Trow = D.T_list(:).';
Tlast = NaN;
iLast = find(det, 1, 'last');
if ~isempty(iLast), Tlast = Trow(iLast); end
end

function [xLo, xHi] = zoom_window(Tlast)
t = Tlast(isfinite(Tlast));
xLo = floor((min(t) - 0.04) / 0.05) * 0.05;
xHi = max(t) + 0.025;
end

function tk = inset_ticks(xLo, xHi)
for step = [0.02 0.05 0.1 0.2 0.5]
    tk = step * (ceil(xLo/step - 1e-9) : floor(xHi/step + 1e-9));
    if numel(tk) <= 3, return; end
end
tk = [xLo xHi];
end

function plot_criteria(D, points, titles, fileBase)
FONT = 'Times New Roman';
critKey = {'ENTSDP', 'ENTCMC', 'ENTPPT'};
critLbl = {'Moment hierarchy', 'CMC', 'PPT witness'};
critCol = [0.000 0.447 0.741; 0.850 0.325 0.098; 0.231 0.549 0.165];
critSty = {'-', '--', '-.'};
yLimits = [1e-4 3];
xPad = 0.08;
insetBox = [0.17 0.15 0.47 0.38];
bcTags = {'OBC', 'PBC'};
nRows = numel(D);
nCols = size(points, 1);
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18.5 1.6 + 4.6*nRows]);
mL = 0.075; mR = 0.012; mB = 0.105; mT = 0.125; gapX = 0.014; gapY = 0.075;
pW = (1 - mL - mR - (nCols - 1)*gapX) / nCols;
pH = (1 - mB - mT - (nRows - 1)*gapY) / nRows;
legendHandles = gobjects(1, numel(critKey));
for ir = 1:nRows
    Trow = D{ir}.T_list(:).';
    for jc = 1:nCols
        ic = case_index(D{ir}, points(jc, :));
        pos = [mL + (jc - 1)*(pW + gapX), mB + (nRows - ir)*(pH + gapY), pW, pH];
        ax = axes(fig, 'Position', pos);
        hold(ax, 'on');
        Y = nan(numel(critKey), numel(Trow));
        Tlast = nan(1, numel(critKey));
        for k = 1:numel(critKey)
            [Y(k, :), Tlast(k)] = detected_series(D{ir}, critKey{k}, ic);
        end
        for k = 1:numel(critKey)
            hk = plot(ax, Trow, Y(k, :), critSty{k}, 'Color', critCol(k, :), 'LineWidth', 1.3);
            if ir == 1 && jc == 1, legendHandles(k) = hk; end
            if isfinite(Tlast(k))
                xline(ax, Tlast(k), ':', 'Color', critCol(k, :), 'LineWidth', 0.9, 'HandleVisibility', 'off');
            end
        end
        if any(isfinite(Tlast))
            xMax = max(Tlast) + xPad;
        else
            xMax = max(Trow);
        end
        set(ax, 'YScale', 'log', 'YLim', yLimits, 'XLim', [min(Trow) xMax], 'FontName', FONT, ...
            'FontSize', 10, 'LineWidth', 0.7, 'TickDir', 'out', 'Box', 'on', 'Layer', 'top');
        grid(ax, 'on');
        xt = 0.5:0.5:3;
        xticks(ax, xt(xt <= xMax));
        if jc > 1, set(ax, 'YTickLabel', []); end
        if ir == 1
            title(ax, titles{jc}, 'Interpreter', 'latex', 'FontSize', 10, 'FontWeight', 'normal');
        end
        if ir == nRows
            xlabel(ax, '$T$', 'Interpreter', 'latex', 'FontSize', 10);
        end
        if jc == 1
            ylabel(ax, [bcTags{ir} ': detection strength'], 'FontName', FONT, 'FontSize', 10);
        end
        if ~any(isfinite(Tlast)), continue; end
        [xLo, xHi] = zoom_window(Tlast);
        inWin = (Trow >= xLo - 1e-12) & (Trow <= xHi + 1e-12);
        v = Y(:, inWin);
        v = v(isfinite(v));
        if isempty(v), continue; end
        yIn = 10.^[floor(log10(min(v))), ceil(log10(max(v)))];
        if yIn(2) <= yIn(1), yIn(2) = 10*yIn(1); end
        fy = [max(yIn(1), yLimits(1)), min(yIn(2), yLimits(2))];
        plot(ax, [xLo xHi xHi xLo xLo], fy([1 1 2 2 1]), '-', 'Color', [0.5 0.5 0.5], ...
            'LineWidth', 0.5, 'HandleVisibility', 'off');
        axIn = axes(fig, 'Position', [pos(1) + insetBox(1)*pos(3), pos(2) + insetBox(2)*pos(4), ...
            insetBox(3)*pos(3), insetBox(4)*pos(4)]);
        hold(axIn, 'on');
        if all(isfinite(Tlast(2:3))) && Tlast(2) > Tlast(3)
            patch(axIn, [Tlast(3) Tlast(2) Tlast(2) Tlast(3)], yIn([1 1 2 2]), [0.88 0.88 0.88], ...
                'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        for k = 1:numel(critKey)
            plot(axIn, Trow, Y(k, :), critSty{k}, 'Color', critCol(k, :), 'LineWidth', 0.9, ...
                'Marker', '.', 'MarkerSize', 7);
            if isfinite(Tlast(k))
                xline(axIn, Tlast(k), ':', 'Color', critCol(k, :), 'LineWidth', 0.8);
            end
        end
        set(axIn, 'YScale', 'log', 'XLim', [xLo xHi], 'YLim', yIn, 'FontName', FONT, 'FontSize', 7, ...
            'LineWidth', 0.5, 'TickDir', 'out', 'TickLength', [0.03 0.03], 'Box', 'on', ...
            'Layer', 'top', 'Color', 'w');
        xticks(axIn, inset_ticks(xLo, xHi));
        dec = log10(yIn(1)):log10(yIn(2));
        if numel(dec) > 3, dec = fliplr(dec(end:-2:1)); end
        yticks(axIn, 10.^dec);
    end
end
lgd = legend(legendHandles, critLbl, 'Orientation', 'horizontal', 'Box', 'off', ...
    'FontName', FONT, 'FontSize', 10);
set(lgd, 'Units', 'normalized');
lp = get(lgd, 'Position');
set(lgd, 'Position', [0.5 - lp(3)/2, 1 - lp(4) - 0.005, lp(3), lp(4)]);
exportgraphics(fig, [fileBase '.pdf'], 'ContentType', 'vector');
end

function plot_gme(D, points, fileBase)
FONT = 'Times New Roman';
critKey = {'ENTSDP', 'ENTCMC', 'ENTPPT'};
mixCol = [0.494 0.184 0.556];
yLim = [1e-4 0.3];
xPad = 0.08;
bcTags = {'OBC', 'PBC'};
nRows = numel(D);
nCols = size(points, 1);
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18.5 1.4 + 4.1*nRows]);
tl = tiledlayout(fig, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
hMix1 = gobjects(1);
for ir = 1:nRows
    Trow = D{ir}.T_list(:).';
    for jc = 1:nCols
        ic = case_index(D{ir}, points(jc, :));
        ax = nexttile(tl);
        hold(ax, 'on');
        [yGME, Tgme] = detected_series(D{ir}, 'FPPTwitStrength', ic);
        Tent = -Inf;
        for k = 1:numel(critKey)
            [~, Tk] = detected_series(D{ir}, critKey{k}, ic);
            Tent = max(Tent, Tk);
        end
        hMix = plot(ax, Trow, yGME, '-', 'Color', mixCol, 'LineWidth', 1.4);
        if ir == 1 && jc == 1, hMix1 = hMix; end
        if isfinite(Tgme)
            xline(ax, Tgme, ':', 'Color', mixCol, 'LineWidth', 1.0, 'HandleVisibility', 'off');
            xMax = Tgme + xPad;
        else
            xMax = max(Trow);
        end
        if isfinite(Tent)
            text(ax, 0.03, 0.045, {'entanglement detected', sprintf('up to {\\itT} = %.2f', Tent)}, ...
                'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
                'FontName', FONT, 'FontSize', 7.5, 'Color', [0.25 0.25 0.25]);
        end
        set(ax, 'YScale', 'log', 'YLim', yLim, 'FontName', FONT, 'FontSize', 10, ...
            'LineWidth', 0.7, 'TickDir', 'out', 'Box', 'on', 'Layer', 'top');
        grid(ax, 'on');
        xlim(ax, [min(Trow) xMax]);
        xt = 0.2:0.2:3;
        xticks(ax, xt(xt <= xMax + 1e-9));
        if jc > 1, set(ax, 'YTickLabel', []); end
        if ir == 1
            title(ax, sprintf('$(\\gamma,h)=(%g,\\,%g)$', points(jc, 1), points(jc, 2)), ...
                'Interpreter', 'latex', 'FontSize', 10, 'FontWeight', 'normal');
        end
        if ir == nRows
            xlabel(ax, '$T$', 'Interpreter', 'latex', 'FontSize', 10);
        end
        if jc == 1
            ylabel(ax, [bcTags{ir} ': detection strength'], 'FontName', FONT, 'FontSize', 10);
        end
    end
end
lgd = legend(hMix1, {'PPT mixer (GME)'}, 'Orientation', 'horizontal', 'Box', 'off', ...
    'FontName', FONT, 'FontSize', 10);
lgd.Layout.Tile = 'north';
exportgraphics(fig, [fileBase '.pdf'], 'ContentType', 'vector');
end
