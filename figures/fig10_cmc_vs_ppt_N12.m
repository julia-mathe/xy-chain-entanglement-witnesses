% FIG10_CMC_VS_PPT_N12  Fig. 10: CMC and PPT witness (alternating cut) in the (gamma, h) plane,
% N = 12, T = 1.47, PBC and OBC.
% Needs CVX. Slow (900 Gibbs states of 12 spins per boundary condition). replot_only = true redraws.

repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));

outdir = fullfile(repo, 'output', 'fig10');
if ~exist(outdir,'dir'), mkdir(outdir); end

N          = 12;
T          = 1.47;
J          = 1.0;
grid_size  = 30;

replot_only = false;
h_max      = 1.5;
gamma_max  = 1.5;

gamma_list = linspace(0, gamma_max, grid_size);
h_list     = linspace(0, h_max,    grid_size);

B_sites    = 2:2:N;

zTol_cmc   = 1e-4;
zTol_ppt   = 1e-6;

data_file = fullfile(outdir, sprintf('cmc_vs_ppt_N%d_T%.3g.mat', N, T));
if replot_only
    load(data_file, 'gamma_list','h_list','T','N','J','omega_pbc','omega_obc', ...
         'ppt_pbc','ppt_obc','B_sites','grid_size');
    fprintf('Loaded saved data: %s\n', data_file);
else
omega_pbc = nan(grid_size, grid_size);
omega_obc = nan(grid_size, grid_size);
ppt_pbc   = nan(grid_size, grid_size);
ppt_obc   = nan(grid_size, grid_size);

ops = local_pauli_operators(N);

fprintf('=== CMC vs PPT phase diagram, N=%d, T=%.3g ===\n', N, T);
fprintf('    Grid: %d x %d = %d points per BC\n', grid_size, grid_size, grid_size^2);

for ibc = 1:2
    is_pbc = (ibc == 1);
    bc_str = 'PBC'; if ~is_pbc, bc_str = 'OBC'; end
    fprintf('\n--- %s ---\n', bc_str);
    tic;

    for ig = 1:grid_size
        gamma = gamma_list(ig);

        for ih = 1:grid_size
            h = h_list(ih);

            H_mat  = xy_hamiltonian(N, gamma, h, is_pbc, J);
            rho    = gibbs_state(H_mat, T);
            [Gam, ~] = covariance_matrix(rho, ops);
            [omega, ~] = cmc_sdp(Gam);

            if is_pbc
                omega_pbc(ih,ig) = real(omega);
            else
                omega_obc(ih,ig) = real(omega);
            end

            rho_PT = partial_transpose(rho, B_sites, N);
            ev_PT  = eig((rho_PT + rho_PT')/2);
            neg_val = min(real(ev_PT));

            if is_pbc
                ppt_pbc(ih,ig) = neg_val;
            else
                ppt_obc(ih,ig) = neg_val;
            end
        end

        fprintf('  %s  gamma=%.2f  (%d/%d, %.0f s)\n', bc_str, gamma, ig, grid_size, toc);
    end
end

save(data_file, ...
    'gamma_list','h_list','T','N','J','omega_pbc','omega_obc','ppt_pbc','ppt_obc', ...
    'B_sites','grid_size','-v7.3');
fprintf('\nData saved.\n');
end

FONT  = 'Times New Roman';
FLAB  = 16;  FTICK = 13;  FTITL = 15;
t_circ  = linspace(0, pi/2, 300);
circ_x  = sin(t_circ);  circ_y = cos(t_circ);

panels = {
    omega_pbc, 'CMC,  PBC',  '$\max(0,1-\omega_{\rm CMC})$', 'cmc', true
    ppt_pbc,   'PPT,  PBC',  '$-\lambda_{\rm min}$',    'ppt', false
    omega_obc, 'CMC,  OBC',  '$\max(0,1-\omega_{\rm CMC})$', 'cmc', true
    ppt_obc,   'PPT,  OBC',  '$-\lambda_{\rm min}$',    'ppt', false
};

fig = figure('Color','w','Units','centimeters','Position',[2 2 26 22]);
tl  = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','loose');

for ip = 1:4
    raw      = real(panels{ip,1});
    ptitle   = panels{ip,2};
    cblabel  = panels{ip,3};
    kind     = panels{ip,4};
    is_cmc   = panels{ip,5};

    ax = nexttile(tl, ip);

    if is_cmc
        detected  = raw < (1 - zTol_cmc);
        plotD     = max(0, 1 - raw);
        cbar_sign = +1;
    else
        detected  = raw < -zTol_ppt;
        plotD     = max(0, -raw);
        cbar_sign = +1;
    end

    maxVal = max(plotD(detected));
    if isempty(maxVal) || maxVal <= 0, maxVal = 1; end

    imagesc(ax, gamma_list, h_list, plotD);
    set(ax,'YDir','normal');
    colormap(ax, parula(256));
    caxis(ax, [0, maxVal]);

    hold(ax,'on');
    rgbBlack = zeros([size(plotD) 3]);
    hBlk     = image(ax, gamma_list, h_list, rgbBlack);
    set(hBlk, 'AlphaData', double(~detected), 'AlphaDataMapping','none');
    uistack(hBlk, 'top');

    plot(ax, circ_x, circ_y, 'w--', 'LineWidth',1.3);
    plot(ax, [0 gamma_max],[1 1], 'w:', 'LineWidth',1.1);
    plot(ax, 0,0,'wo','MarkerSize',7,'MarkerFaceColor','w');
    plot(ax, 1,1,'ws','MarkerSize',7,'MarkerFaceColor','w');

    cb = colorbar(ax);
    cb.FontName = FONT; cb.FontSize = FTICK;
    ylabel(cb, cblabel, 'Interpreter','latex','FontSize',FLAB);

    set(ax,'FontName',FONT,'FontSize',FTICK,'TickDir','out','Box','on', ...
        'XLim',[0 gamma_max],'YLim',[0 h_max]);
    xlabel(ax,'$\gamma$','Interpreter','latex','FontSize',FLAB);
    ylabel(ax,'$h/J$',   'Interpreter','latex','FontSize',FLAB);
    title(ax, ptitle, 'FontName',FONT,'FontSize',FTITL,'FontWeight','normal');

    text(ax,0.03,0.97,'Black: not detected','Units','normalized', ...
        'VerticalAlignment','top','Color','w','FontName',FONT,'FontSize',11);
end

sgtitle(fig, sprintf('CMC vs PPT,  $N=%d$,  $T/J=%.3g$', N, T), ...
    'Interpreter','latex','FontSize',FTITL);

exportgraphics(fig, fullfile(outdir, sprintf('cmc_vs_ppt_N%d_T%.3g.pdf',N,T)), ...
    'ContentType','vector');
exportgraphics(fig, fullfile(outdir, sprintf('cmc_vs_ppt_N%d_T%.3g.png',N,T)), ...
    'Resolution',300);
fprintf('Figure saved.\n');
