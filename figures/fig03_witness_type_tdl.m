% FIG03_WITNESS_TYPE_TDL  Fig. 3: omega_CMC and optimal witness type of the infinite chain, T = 1e-3.
% Set run_obc = true for the window of 80 spins in an open chain of 100 (Sec. 4.1; needs CVX, slow);
% use_saved_data = true redraws from output/fig03.

repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));

outdir = fullfile(repo, 'output', 'fig03');
if ~exist(outdir, 'dir'), mkdir(outdir); end

T          = 1e-3;
grid_size  = 25;
gamma_max  = 1.5;
h_max      = 1.5;
gamma_list = linspace(0, gamma_max, grid_size);
h_list     = linspace(0, h_max,    grid_size);

run_obc    = true;
N_big      = 100;
N_sub      = 80;
solverName = 'mosek';

use_saved_data = false;

save_path = fullfile(outdir, sprintf('large_N_map_T%.0e.mat', T));
if use_saved_data
    load(save_path, 'gamma_list','h_list','T','omega_tdl','omega_sub', ...
         's_tdl','q_tdl','N_big','N_sub','grid_size');
    fprintf('Loaded saved data: %s\n', save_path);
else
fprintf('=== Panel (2): thermodynamic limit (%dx%d grid, T=%.0e) ===\n', ...
    grid_size, grid_size, T);
omega_tdl = zeros(grid_size, grid_size);
s_tdl     = zeros(grid_size, grid_size, 3);
q_tdl     = zeros(grid_size, grid_size, 3);
tic;
for ig = 1:grid_size
    for ih = 1:grid_size
        [omega_tdl(ih,ig), s_ij, q_ij] = cmc_tdl(gamma_list(ig), h_list(ih), T);
        s_tdl(ih,ig,:) = s_ij;
        q_tdl(ih,ig,:) = q_ij;
    end
end
fprintf('  done (%.0f s)\n', toc);

omega_sub = NaN(grid_size, grid_size);

if run_obc
    center     = floor(N_big/2);
    I_win      = (center - floor(N_sub/2)) : (center + ceil(N_sub/2) - 1);
    I_win      = I_win(I_win >= 1 & I_win <= N_big);
    N_sub_act  = length(I_win);
    first_site = I_win(1);

    fprintf('=== Panel (1): OBC subchain N_sub=%d in N_big=%d (%dx%d grid) ===\n', ...
        N_sub_act, N_big, grid_size, grid_size);
    fprintf('  Solver: %s.  SDP size: %dx%d.\n', solverName, 3*N_sub_act, 3*N_sub_act);
    fprintf('  Estimated time: %.0f-%.0f min (MOSEK) or %.0f-%.0f min (SDPT3).\n', ...
        grid_size^2*1/60, grid_size^2*3/60, grid_size^2*8/60, grid_size^2*15/60);

    tic;
    for ig = 1:grid_size
        for ih = 1:grid_size
            G_big = majorana_G_open(N_big, gamma_list(ig), h_list(ih), T);
            Gam   = covariance_window_open(G_big, first_site, N_sub_act);
            omega_sub(ih,ig) = cmc_sdp(Gam);
        end
        fprintf('  gamma = %.2f  (%d/%d, %.0f s)\n', gamma_list(ig), ig, grid_size, toc);
    end
    fprintf('  OBC panel done (%.0f s)\n', toc);
end

save(save_path, 'gamma_list','h_list','T','omega_tdl','omega_sub', ...
     's_tdl','q_tdl','N_big','N_sub','grid_size', '-v7.3');
fprintf('Data saved: %s\n', save_path);
end

FONT      = 'Times New Roman';
FLAB      = 16;
FTICK     = 13;
FTITL     = 15;
detectTol = 1e-4;
t_circ    = linspace(0, pi/2, 300);

circ_x = sin(t_circ);   circ_y = cos(t_circ);
ising_x = [0 gamma_max]; ising_y = [1 1];

if run_obc
    fprintf('OBC window: %d of %d grid points not detected\n', nnz(omega_sub >= 1), numel(omega_sub));
end

FLAB      = 18;
FTICK     = 18;
FTITL     = 18;
cat_colors = [ 31  90 166;
              209  73  91;
              143 192 232;
              242 177  52;
              122 158  59;
              148  52 198;
              217 217 217;
              110 110 110] / 255;

cat_names  = {'$\mathcal{S}_x(\pi),\mathcal{S}_y(\pi),\mathcal{S}_z(0)$ (XX)', ...
              '$\mathcal{S}_x(\pi),\mathcal{S}_y(0),\mathcal{S}_z(\pi)$ (Ising)', ...
              '$\mathcal{S}_x(\pi),\mathcal{S}_z(0)$', ...
              '$\mathcal{S}_x(\pi),\mathcal{S}_z(\pi)$', ...
              '$\mathcal{S}_y(0),\mathcal{S}_z(\pi)$', ...
              'incommensurate minimum', ...
              '$\omega_{\rm CMC}>0.999$ (undetected)', ...
              'other commensurate witness'};

eps_q   = 0.15;
eps_tie = 0.01;

catmap     = zeros(grid_size, grid_size);
tied_cats  = cell(grid_size, grid_size);

for ig = 1:grid_size
    for ih = 1:grid_size
        q_ij = squeeze(q_tdl(ih,ig,:))';
        s_ij = squeeze(s_tdl(ih,ig,:))';
        is0  = abs(q_ij)      < eps_q;
        ispi = abs(q_ij - pi) < eps_q;

        if omega_tdl(ih,ig) > 0.999
            catmap(ih,ig) = 7; continue;
        end

        sx = s_ij(1); sy = s_ij(2); sz = s_ij(3);
        cands    = [(sx+sy+sz)/2, sx+sy, sy+sz, sz+sx];
        cand_min = min(cands);
        tied     = abs(cands - cand_min) < eps_tie * max(abs(cand_min), 1e-10);
        n_tied   = sum(tied);

        comps = {[1 2 3], [1 2], [2 3], [1 3]};
        cand_to_cat = 8*ones(1,4);
        for c = 1:4
            a = comps{c};
            if ~all(is0(a) | ispi(a))
                cand_to_cat(c) = 6;
            elseif c == 1 && ispi(1) && ispi(2) && is0(3)
                cand_to_cat(c) = 1;
            elseif c == 1 && ispi(1) && is0(2) && ispi(3)
                cand_to_cat(c) = 2;
            elseif c == 3 && is0(2) && ispi(3)
                cand_to_cat(c) = 5;
            elseif c == 4 && ispi(1) && is0(3)
                cand_to_cat(c) = 3;
            elseif c == 4 && ispi(1) && ispi(3)
                cand_to_cat(c) = 4;
            end
        end

        [~, best] = min(cands);
        catmap(ih,ig) = cand_to_cat(best);

        if n_tied > 1
            tc = unique(cand_to_cat(tied));
            tc = tc(tc <= 6);
            if numel(tc) > 1
                tied_cats{ih,ig} = tc;
            end
        end
    end
end

fig3   = figure('Color','w','Units','centimeters','Position',[2 2 38 16]);
ax_bot = 0.12;
ax_h   = 0.8;
ax_w   = 0.3;
ax_gap = 0.1;
ax_L0  = 0.06;

axL = axes('Parent',fig3,'Position',[ax_L0,             ax_bot, ax_w, ax_h]);
axR = axes('Parent',fig3,'Position',[ax_L0+ax_w+ax_gap, ax_bot, ax_w, ax_h]);

axes(axL);
imagesc(axL, gamma_list, h_list, omega_tdl);
set(axL, 'YDir','normal');
colormap(axL, viridis_colormap(256));
caxis(axL, [0.1 1]);
cb = colorbar(axL);
ylabel(cb, '$\omega_{\rm CMC}$','Interpreter','latex','FontSize',FLAB);
cb.FontName = FONT;  cb.FontSize = FTICK;
hold(axL,'on');
[cm,ch_c] = contour(axL, gamma_list, h_list, omega_tdl, [0.25 0.50 0.75 0.95], ...
    'LineColor','w', 'LineWidth',1.0);
clabel(cm, ch_c, 'FontSize',10, 'Color','w', 'FontName',FONT);
mask_L = omega_tdl >= (1 - detectTol);
h_blk  = image(axL, gamma_list, h_list, zeros([size(mask_L) 3]));
set(h_blk, 'AlphaData', double(mask_L), 'AlphaDataMapping','none');
uistack(h_blk,'top');
plot(axL, circ_x,  circ_y,  'w--', 'LineWidth',1.4);
plot(axL, ising_x, ising_y, 'w:',  'LineWidth',1.2);
plot(axL, 0,0,'wo','MarkerSize',8,'MarkerFaceColor','w');
plot(axL, 1,1,'ws','MarkerSize',8,'MarkerFaceColor','w');
set(axL,'FontName',FONT,'FontSize',FTICK,'TickDir','out','Box','on', ...
    'XLim',[0 gamma_max],'YLim',[0 h_max]);
xlabel(axL,'$\gamma$','Interpreter','latex','FontSize',FLAB);
ylabel(axL,'$h/J$','Interpreter','latex','FontSize',FLAB);
T_exp   = floor(log10(T));
T_coeff = T / 10^T_exp;
if abs(T_coeff - 1) < 1e-9
    T_str = sprintf('10^{%d}', T_exp);
else
    T_str = sprintf('%.1f\\times 10^{%d}', T_coeff, T_exp);
end
title(axL, sprintf('(a) $\\omega_{\\rm CMC}$,\\quad$N\\to\\infty$,\\quad$T/J=%s$', T_str), ...
    'Interpreter','latex','FontSize',FTITL,'FontWeight','normal');

axes(axR);

pic = zeros([size(catmap) 3]);
for ch = 1:3
    col_channel = cat_colors(:,ch);
    pic(:,:,ch) = col_channel(catmap);
end

up  = 8;
Nh  = grid_size; Ng = grid_size;
pic_up = zeros(Nh*up, Ng*up, 3);
for ch = 1:3
    col_channel = cat_colors(:,ch);
    tmp = col_channel(catmap);
    pic_up(:,:,ch) = repelem(tmp, up, up);
end
stripe_w = max(1, floor(up/2));
for ig = 1:Ng
    for ih = 1:Nh
        tc = tied_cats{ih,ig};
        if isempty(tc), continue; end
        n_tc = numel(tc);
        r0 = (ih-1)*up + 1;  r1 = ih*up;
        c0 = (ig-1)*up + 1;  c1 = ig*up;
        for rr = r0:r1
            for cc = c0:c1
                diag_idx = mod(floor((rr-r0 + cc-c0) / stripe_w), n_tc) + 1;
                cat_here = tc(diag_idx);
                for ch = 1:3
                    pic_up(rr,cc,ch) = cat_colors(cat_here,ch);
                end
            end
        end
    end
end

image(axR, gamma_list, h_list, pic_up);
set(axR, 'YDir','normal');
hold(axR,'on');
plot(axR, circ_x,  circ_y,  'k--','LineWidth',1.4);
plot(axR, ising_x, ising_y, 'k:', 'LineWidth',1.2);
plot(axR, 0,0,'ko','MarkerSize',8,'MarkerFaceColor','k');
plot(axR, 1,1,'ks','MarkerSize',8,'MarkerFaceColor','k');
set(axR,'FontName',FONT,'FontSize',FTICK,'TickDir','out','Box','on', ...
    'XLim',[0 gamma_max],'YLim',[0 h_max]);
xlabel(axR,'$\gamma$','Interpreter','latex','FontSize',FLAB);
ylabel(axR,'$h/J$','Interpreter','latex','FontSize',FLAB);
title(axR,'(b) Optimal witness type, $N\to\infty$', ...
    'Interpreter','latex','FontSize',FTITL,'FontWeight','normal');

short_names = {'$\mathcal{S}_x(\pi),\mathcal{S}_y(\pi),\mathcal{S}_z(0)$', ...
               '$\mathcal{S}_x(\pi),\mathcal{S}_y(0),\mathcal{S}_z(\pi)$', ...
               '$\mathcal{S}_x(\pi),\mathcal{S}_z(0)$', ...
               '$\mathcal{S}_x(\pi),\mathcal{S}_z(\pi)$', ...
               '$\mathcal{S}_y(0),\mathcal{S}_z(\pi)$', ...
               'incommensurate', ...
               '$\omega_{\rm CMC}>0.999$', ...
               'other'};

present_cats = unique(catmap(:))';
for ig2 = 1:grid_size
    for ih2 = 1:grid_size
        present_cats = union(present_cats, tied_cats{ih2,ig2});
        present_cats = present_cats(:)';
    end
end
present_primary = [present_cats(present_cats ~= 7 & present_cats ~= 8), ...
                   present_cats(present_cats == 8), present_cats(present_cats == 7)];

tied_combo_list = {};
for ig2 = 1:Ng
    for ih2 = 1:Nh
        tc = tied_cats{ih2,ig2};
        if isempty(tc), continue; end
        key = sprintf('%d+',tc);
        if ~any(strcmp(tied_combo_list, key))
            tied_combo_list{end+1} = key;
        end
    end
end
tied_combo_list = sort(tied_combo_list);

n_prim  = numel(present_primary);
n_combo = numel(tied_combo_list);
n_total = n_prim + n_combo;

axDummy = axes('Parent',fig3,'Visible','off','Position',[0.001 0.001 0.001 0.001]);
hold(axDummy,'on');
h_proxy   = gobjects(n_total,1);
lgd_names = cell(n_total,1);

for k = 1:n_prim
    c = present_primary(k);
    h_proxy(k) = patch(axDummy, NaN, NaN, cat_colors(c,:), ...
        'EdgeColor','k','LineWidth',0.5);
    lgd_names{k} = short_names{c};
end

lgd = legend(axDummy, h_proxy(1:n_prim), lgd_names(1:n_prim), ...
    'Interpreter','latex','FontName',FONT,'FontSize',FTITL, ...
    'Box','off','NumColumns',1);
lgd.Units    = 'normalized';
lgd.Position = [0.72, 0.06, 0.27, 0.88];

exportgraphics(fig3, fullfile(outdir,'witness_type_TDL.pdf'), 'ContentType','vector');
exportgraphics(fig3, fullfile(outdir,'witness_type_TDL.png'), 'Resolution',300);
fprintf('Saved: witness_type_TDL\n');
fprintf('\nAll figures in: %s\n', outdir);
