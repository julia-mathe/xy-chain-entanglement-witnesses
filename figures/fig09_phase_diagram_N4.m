% FIG09_PHASE_DIAGRAM_N4  Fig. 9: moment hierarchy, CMC and PPT witness in the (gamma, h) plane,
% N = 4, T = 1.47. Run once with BC = "open" and once with BC = "periodic".
% Needs CVX. (MOSEK recommended). Slow. replot_only = true redraws from output/fig09.

repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));

Nq    = 4;
BC    = "periodic";
J     = 1.0;
T     = 1.47;
beta  = 1 / T;
B     = 2:2:Nq;

discretization = 30;
h_list     = linspace(0, 1.5, discretization);
gamma_list = linspace(0, 1.5, discretization);

kDual      = 2;
solverName = 'mosek';
zeroTol    = 1e-8;

replot_only = false;

paramString = sprintf('N%d_T%.4g_J%.4g_%s_k%d', Nq, T, J, char(BC), kDual);
paramString = strrep(strrep(paramString,'.','p'),'-','m');

outputDir = fullfile(repo, 'output', 'fig09', paramString);
figDir    = fullfile(outputDir, 'figures');
if ~exist(figDir,'dir'), mkdir(figDir); end
dataFile  = fullfile(outputDir, ['results_' paramString '.mat']);

fprintf('\n============================================================\n');
fprintf(' PHASE DIAGRAM SCAN: moment hierarchy / CMC / PPT\n');
fprintf(' N=%d, T=%.4g, J=%.4g, BC=%s, k=%d, solver=%s\n', ...
    Nq, T, J, char(BC), kDual, solverName);
fprintf(' Grid: %d x %d = %d points\n', ...
    discretization, discretization, discretization^2);
fprintf(' Output: %s\n', outputDir);
fprintf('============================================================\n\n');

if replot_only
    load(dataFile, 'settings','ENTSDP','ENTCMC','ENTPPT', ...
         'rawSDP','rawCMC','rawPPT');
    h_list = settings.h_list;  gamma_list = settings.gamma_list;
    Nh = numel(h_list);  Ng = numel(gamma_list);
    fprintf('Loaded saved data: %s\n', dataFile);
else
Nh = numel(h_list);
Ng = numel(gamma_list);

ENTSDP    = nan(Nh, Ng);
ENTCMC    = nan(Nh, Ng);
ENTPPT    = nan(Nh, Ng);

rawSDP    = nan(Nh, Ng);
rawCMC    = nan(Nh, Ng);
rawPPT    = nan(Nh, Ng);

ops   = local_pauli_operators(Nq);

tic;
for u = 1:Nh
    hval = h_list(u);
    fprintf('h = %.4f  (%d/%d, %.0f s elapsed)\n', hval, u, Nh, toc);

    for d = 1:Ng
        gamma = gamma_list(d);

        rho = gibbs_state(xy_hamiltonian(Nq, gamma, hval, BC == "periodic", J), T);
        v = moment_hierarchy_value(rho, Nq, kDual, solverName);
        rawSDP(u,d) = v;
        ENTSDP(u,d) = max(0, -v);

        Gamma        = covariance_matrix(rho, ops);
        [omega, ~]   = cmc_sdp(Gamma);
        rawCMC(u,d)  = real(omega);
        ENTCMC(u,d)  = max(0, 1 - real(omega));

        neg_val = ppt_negative_part(rho, B, Nq);
        rawPPT(u,d)  = real(neg_val);
        ENTPPT(u,d)  = max(0, -real(neg_val));
    end
end
fprintf('\nScan complete (%.0f s total)\n\n', toc);

settings = struct('Nq',Nq,'BC',char(BC),'J',J,'T',T,'beta',beta,'B',B, ...
    'discretization',discretization,'h_list',h_list,'gamma_list',gamma_list, ...
    'kDual',kDual,'solverName',solverName,'zeroTol',zeroTol);
save(dataFile,'settings','ENTSDP','ENTCMC','ENTPPT', ...
    'rawSDP','rawCMC','rawPPT','-v7.3');
fprintf('Data saved: %s\n', dataFile);
end

FONT      = 'Times New Roman';
FLAB      = 18;
FTICK     = 18;
FTITL     = 18;
FCBAR     = 18;
t_circ    = linspace(0, pi/2, 300);
plotMap   = parula(256);

fmt = @(ax) set(ax, 'FontName',FONT, 'FontSize',FTICK, ...
    'LineWidth',1.2, 'TickDir','out', 'Box','on', 'Layer','top');

data_list  = {ENTSDP,   ENTCMC,   ENTPPT};
name_list  = {'Moment hierarchy', 'CMC', 'PPT witness'};
cbar_list  = {'$\max(0,-v_{\mathrm{SDP}})$', ...
              '$\max(0,1-\omega_{\mathrm{CMC}})$', ...
              '$\max(0,-v_{\mathrm{PPT}})$'};
file_list  = {'SDP','CMC','PPT'};

fig_comb = figure('Color','w','Units','centimeters','Position',[2 2 38 13]);
tl = tiledlayout(fig_comb, 1, 3, 'TileSpacing','compact', 'Padding','loose');

for iq = 1:3
    Q      = real(data_list{iq});
    finQ   = Q(isfinite(Q));
    maxQ   = max([finQ(:); 0]);
    cLim   = [0, max(maxQ, eps)];

    ax = nexttile(tl, iq);

    plotData  = Q;
    alphaMask = double(isfinite(plotData) & plotData >= zeroTol);

    hImg = imagesc(ax, gamma_list, h_list, plotData);
    hold(ax, 'on');
    rgbBlack = zeros([size(plotData), 3]);
    hBlack   = image(ax, gamma_list, h_list, rgbBlack);
    set(hBlack, 'AlphaData', 1 - alphaMask, 'AlphaDataMapping','none');
    uistack(hBlack, 'top');

    plot(ax, sin(t_circ), cos(t_circ), 'w--', 'LineWidth',1.4);
    plot(ax, [0 max(gamma_list)], [1 1], 'w:', 'LineWidth',1.2);
    plot(ax, 0, 0, 'wo', 'MarkerSize',8, 'MarkerFaceColor','w');
    plot(ax, 1, 1, 'ws', 'MarkerSize',8, 'MarkerFaceColor','w');

    set(ax, 'YDir','normal');
    colormap(ax, plotMap);
    caxis(ax, cLim);
    fmt(ax);

    xlabel(ax, '$\gamma$', 'Interpreter','latex','FontSize',FLAB);
    ylabel(ax, '$h/J$',    'Interpreter','latex','FontSize',FLAB);
    title(ax, sprintf('%s', name_list{iq}), ...
        'Interpreter','latex','FontSize',FTITL,'FontWeight','normal');

    cb = colorbar(ax);
    cb.FontName = FONT;  cb.FontSize = FCBAR;
    cb.Label.String      = cbar_list{iq};
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize    = FCBAR;

    axis(ax,'tight'); axis(ax,'square');
end

exportgraphics(fig_comb, ...
    fullfile(figDir,['phase_diagram_combined_' paramString '.pdf']), ...
    'ContentType','vector');
exportgraphics(fig_comb, ...
    fullfile(figDir,['phase_diagram_combined_' paramString '.png']), ...
    'Resolution',400);
fprintf('Saved: phase_diagram_combined\n');

fprintf('\n============================================================\n');
fprintf(' SCAN COMPLETE\n');
fprintf('============================================================\n');
nGrid = Nh * Ng;
for iq = 1:3
    Q = real(data_list{iq});
    nDet = nnz(isfinite(Q) & Q >= zeroTol);
    fprintf('%s: %d / %d grid points detected (%.1f%%)\n', ...
        name_list{iq}, nDet, nGrid, 100*nDet/nGrid);
end
fprintf('All output in: %s\n', outputDir);
fprintf('============================================================\n');
