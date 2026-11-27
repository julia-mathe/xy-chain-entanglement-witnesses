% FIG07_14_QUENCH  Figs. 7 and 14: CMC after a quench from a thermal state, infinite Ising chain.
% Needs CVX. The map of Fig. 14(a) takes a few minutes.

repo = fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(repo,'functions'));
outdir = fullfile(repo, 'output', 'fig07_14'); if ~exist(outdir,'dir'), mkdir(outdir); end

%% SETTINGS
J = 1; gamma_i = 1; h_i = 0.3; gamma_f = 1; h_f = 2; T_i = 1;
t_list  = 0:0.1:10;
nk = 2048;  nq = 256;
Ti_scan = linspace(0.5,2,26);
do_appendix = true;
h_grid  = linspace(0.05,3,25);
L_conv  = 4:4:20;  N_conv = [8 12 16 20 24];
do_ti_sdp = true;

FONT='Times New Roman'; FL=13; FT=11; grey=[0.92 0.92 0.92];
c_red=[209 73 91]/255; c_or=[242 150 40]/255; c_teal=[42 157 143]/255; c_blue=[31 90 166]/255;
c_green=[122 158 59]/255; c_dark=[0.35 0.35 0.35];

ff_selftest(); tdl_selftest();

%% STATES: initial, steady state (GGE), Gibbs state of H_f with the same energy
S_i = tdl_setup(J,gamma_i,h_i,nk);  S_f = tdl_setup(J,gamma_f,h_f,nk);
mu0 = tdl_mu_gibbs(S_i.A,T_i);  mu_gge = tdl_mu_gge(mu0,S_f.A);
T_eff = tdl_teff(S_f.A,tdl_energy(S_f.A,mu0));
[w0,s0] = theorem1_from_corr(tdl_corr_adaptive(mu0,S_i.k));
[wG,sG] = theorem1_from_corr(tdl_corr_adaptive(mu_gge,S_f.k));
[wT,sT] = theorem1_from_corr(tdl_corr_adaptive(tdl_mu_gibbs(S_f.A,T_eff),S_f.k));

%% DYNAMICS
vmax = max(abs(diff(S_f.eps)))/(2*pi/nk);  R = ceil(2*vmax*t_list(end)) + 20;
nt = numel(t_list); w = nan(1,nt); pw = w; b3 = w; tic;
for it = 1:nt
    Bm = tdl_blocks(tdl_mu_evolve(mu0,S_f.A,t_list(it)),S_f.k,R);
    [Gd,mz] = tdl_corr(Bm,R);
    if do_ti_sdp, w(it) = cmc_ti_sdp(Gd,nq); end
    if isnan(w(it)), w(it) = theorem1_from_corr(Gd); end
    if do_appendix, pw(it) = pairwise_min_tdl(Gd,mz,12); b3(it) = block_ppt_min(tdl_window(Bm,3),3); end
    if mod(it,25)==0, fprintf('  t = %.1f (%.0f s)\n',t_list(it),toc); end
end

%% INITIAL-TEMPERATURE WINDOW
RT = nan(numel(Ti_scan),3);
for n = 1:numel(Ti_scan)
    m0 = tdl_mu_gibbs(S_i.A,Ti_scan(n)); Te = tdl_teff(S_f.A,tdl_energy(S_f.A,m0));
    RT(n,:) = [theorem1_from_corr(tdl_corr_adaptive(m0,S_i.k)), ...
               theorem1_from_corr(tdl_corr_adaptive(tdl_mu_gge(m0,S_f.A),S_f.k)), ...
               theorem1_from_corr(tdl_corr_adaptive(tdl_mu_gibbs(S_f.A,Te),S_f.k))];
end
T_ent = crossing(Ti_scan,RT(:,1)');  T_star = crossing(Ti_scan,RT(:,2)');

%% MAIN FIGURE
fig = figure('Color','w','Units','centimeters','Position',[2 2 30 10]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); yl = [min(w)-0.03, max([wT w])+0.05];
shade(ax,[0 t_list(end)],[1 yl(2)],grey); yline(ax,1,'-k','HandleVisibility','off');
yline(ax,wT,':','Color',c_dark,'LineWidth',1.5,'DisplayName','Gibbs, same energy');
yline(ax,wG,'--','Color',c_green,'LineWidth',1.3,'DisplayName','steady state (GGE)');
scatter(ax,t_list,w,16,'filled','MarkerFaceColor',c_red,'DisplayName','$\omega_{\rm CMC}(t)$');
fmt(ax,FONT,FT,[0 t_list(end)],yl);
xlabel(ax,'$tJ$','Interpreter','latex','FontSize',FL); ylabel(ax,'$\omega_{\rm CMC}$','Interpreter','latex','FontSize',FL);
legend(ax,'Interpreter','latex','Box','off','FontSize',FT-1,'NumColumns',3,'Location','north');
title(ax,sprintf('(a) $h_i=%.2gJ\\to h_f=%.2gJ$, $T_i=%.2gJ$',h_i,h_f,T_i),'Interpreter','latex','FontSize',FL,'FontWeight','normal');
ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); yl = [min(RT(:))-0.04, max(RT(:))+0.06];
shade(ax,Ti_scan([1 end]),[1 yl(2)],grey); yline(ax,1,'-k','HandleVisibility','off');
if ~isnan(T_ent) && ~isnan(T_star)
    fill(ax,[T_ent T_star T_star T_ent],[yl(1) yl(1) yl(2) yl(2)],c_red,'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');
end
scatter(ax,Ti_scan,RT(:,1),20,'o','MarkerEdgeColor',c_dark,'DisplayName','initial state');
scatter(ax,Ti_scan,RT(:,2),20,'filled','MarkerFaceColor',c_green,'DisplayName','steady state (GGE)');
scatter(ax,Ti_scan,RT(:,3),20,'filled','MarkerFaceColor',c_dark,'DisplayName','Gibbs, same energy');
fmt(ax,FONT,FT,Ti_scan([1 end]),yl);
xlabel(ax,'$T_i/J$','Interpreter','latex','FontSize',FL); ylabel(ax,'$\omega_{\rm CMC}$','Interpreter','latex','FontSize',FL);
legend(ax,'Interpreter','latex','Box','off','FontSize',FT-1,'Location','southeast');
title(ax,'(b) initial-temperature window','Interpreter','latex','FontSize',FL,'FontWeight','normal');
save_fig(fig,outdir,'quench_main');

%% APPENDIX FIGURE
if do_appendix
    ng = numel(h_grid); WI = nan(ng,1); WGm = nan(ng);
    for a = 1:ng
        Sa = tdl_setup(J,gamma_i,h_grid(a),nk); m0 = tdl_mu_gibbs(Sa.A,T_i);
        WI(a) = theorem1_from_corr(tdl_corr_adaptive(m0,Sa.k));
        for b = 1:ng
            Sb = tdl_setup(J,gamma_f,h_grid(b),nk);
            WGm(a,b) = theorem1_from_corr(tdl_corr_adaptive(tdl_mu_gge(m0,Sb.A),Sb.k));
        end
    end
    BmG = tdl_blocks(mu_gge,S_f.k,max(L_conv));
    wL = arrayfun(@(L) cmc_sdp(spin_covariance_ff(tdl_window(BmG,L),L)), L_conv);
    wN = arrayfun(@(N) cmc_sdp(spin_covariance_ff(majorana_gge(majorana_M_gibbs( ...
             majorana_A(N,J,gamma_i,h_i),T_i),majorana_A(N,J,gamma_f,h_f)),N)), N_conv);

    fig = figure('Color','w','Units','centimeters','Position',[2 2 36 10]);
    tl = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
    ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
    cats = 4*ones(ng); before = WI < 1; after = WGm < 1;
    cats(before & after) = 1; cats(~before & after) = 2; cats(before & ~after) = 3;
    ccol = [143 192 232; 209 73 91; 242 150 40; 235 235 235]/255; pic = zeros(ng,ng,3);
    for ch = 1:3, cc = ccol(:,ch); pic(:,:,ch) = cc(cats'); end
    image(ax,h_grid,h_grid,pic); set(ax,'YDir','normal');
    xline(ax,1,':k'); yline(ax,1,':k'); plot(ax,h_i,h_f,'kp','MarkerFaceColor','k','MarkerSize',10);
    hp = gobjects(4,1); for c = 1:4, hp(c) = patch(ax,NaN,NaN,ccol(c,:),'EdgeColor','k'); end
    legend(ax,hp,{'detected before and after','generated by the quench','lost','never detected'},'FontSize',FT-2,'Location','northeast');
    fmt(ax,FONT,FT,h_grid([1 end]),h_grid([1 end]));
    xlabel(ax,'$h_i/J$','Interpreter','latex','FontSize',FL); ylabel(ax,'$h_f/J$','Interpreter','latex','FontSize',FL);
    title(ax,sprintf('(a) which quench? ($T_i=%.2gJ$, steady state)',T_i),'Interpreter','latex','FontSize',FL,'FontWeight','normal');
    ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); yl = [1.15*min([pw b3]), 1.3*max([pw b3 1e-3])];
    shade(ax,[0 t_list(end)],[0 yl(2)],grey); yline(ax,0,'-k','HandleVisibility','off');
    scatter(ax,t_list,pw,14,'filled','MarkerFaceColor',c_or,'DisplayName','all pairs');
    scatter(ax,t_list,b3,14,'filled','MarkerFaceColor',c_teal,'DisplayName','contiguous 3-spin blocks');
    fmt(ax,FONT,FT,[0 t_list(end)],yl);
    xlabel(ax,'$tJ$','Interpreter','latex','FontSize',FL); ylabel(ax,'$\lambda_{\min}(\rho^{T})$','Interpreter','latex','FontSize',FL);
    legend(ax,'Interpreter','latex','Box','off','FontSize',FT-1,'Location','southeast');
    title(ax,'(b) entanglement of reduced states','Interpreter','latex','FontSize',FL,'FontWeight','normal');
    ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
    yline(ax,wG,'--','Color',c_green,'LineWidth',1.2,'DisplayName','infinite chain (Theorem 1)');
    scatter(ax,L_conv,wL,30,'o','MarkerEdgeColor',c_blue,'LineWidth',1,'DisplayName','window of $L$ spins');
    scatter(ax,N_conv,wN,30,'filled','MarkerFaceColor',c_red,'DisplayName','open chain of $N$ spins');
    fmt(ax,FONT,FT,[0 max([L_conv N_conv])+2],[min([wL wN wG])-0.01, max([wL wN])+0.01]);
    xlabel(ax,'$L$ or $N$','Interpreter','latex','FontSize',FL); ylabel(ax,'$\omega_{\rm CMC}$ (steady state)','Interpreter','latex','FontSize',FL);
    legend(ax,'Interpreter','latex','Box','off','FontSize',FT-1);
    title(ax,'(c) convergence','Interpreter','latex','FontSize',FL,'FontWeight','normal');
    save_fig(fig,outdir,'quench_appendix');
end

%% SUMMARY (numbers quoted in the text)
fprintf('\nomega_CMC: initial %.3f | GGE %.3f | Gibbs(T_eff = %.2fJ) %.3f | first detection t = %.1f/J\n', ...
    w0,wG,T_eff,wT,t_list(find(w<1,1)));
fprintf('window: T_ent(H_i) = %.2fJ < T_i < T_i* = %.2fJ\n',T_ent,T_star);
fprintf('structure factors s_x,s_y,s_z: initial %.2f %.2f %.2f | GGE %.2f %.2f %.2f | Gibbs %.2f %.2f %.2f\n',s0,sG,sT);
if do_appendix, fprintf('min pair eigenvalue %+.4f | 3-blocks entangled at %d of %d times\n',min(pw),sum(b3<0),nt); end

%% ======================= LOCAL FUNCTIONS =======================
function S = tdl_setup(J,g,h,nk)
    S.k = 2*pi*(0:nk-1)/nk;  B0 = [0 2*h; -2*h 0];  B1 = [0 -J*(1-g); J*(1+g) 0];
    S.A = zeros(2,2,nk); S.eps = zeros(1,nk);
    for n = 1:nk
        S.A(:,:,n) = B0 + B1*exp(1i*S.k(n)) - B1.'*exp(-1i*S.k(n));
        S.eps(n) = max(real(eig(1i*S.A(:,:,n))));
    end
end

function mu = tdl_mu_gibbs(A,T)
    mu = zeros(size(A));
    for n = 1:size(A,3)
        X = 1i*A(:,:,n); X = (X+X')/2; [U,L] = eig(X); lam = real(diag(L));
        if T == 0, f = sign(lam); else, f = tanh(lam/(2*T)); end
        mu(:,:,n) = 1i*(U*diag(f)*U');
    end
end

function mu = tdl_mu_evolve(mu0,A,t)
    mu = zeros(size(mu0));
    for n = 1:size(A,3), O = expm(A(:,:,n)*t); mu(:,:,n) = O*mu0(:,:,n)*O'; end
end

function mu = tdl_mu_gge(mu0,A)
    mu = zeros(size(mu0));
    for n = 1:size(A,3)
        X = 1i*A(:,:,n); X = (X+X')/2; [U,~] = eig(X);
        mu(:,:,n) = U*diag(diag(U'*mu0(:,:,n)*U))*U';
    end
end

function e = tdl_energy(A,mu)
    e = 0; for n = 1:size(A,3), e = e - 0.25*real(trace(A(:,:,n)*mu(:,:,n))); end
    e = e/size(A,3);
end

function T = tdl_teff(A,E)
    lo = 1e-3; hi = 1e3;
    for it = 1:80
        mid = sqrt(lo*hi);
        if tdl_energy(A,tdl_mu_gibbs(A,mid)) < E, lo = mid; else, hi = mid; end
    end
    T = sqrt(lo*hi);
end

function Bm = tdl_blocks(mu,k,R)
    nk = numel(k); E = exp(-1i*k(:)*(-R:R));
    Bm = real(reshape(reshape(mu,4,nk)*E/nk,2,2,2*R+1));
end

function M = tdl_window(Bm,L)
    R = (size(Bm,3)-1)/2; M = zeros(2*L);
    for j = 1:L, for l = 1:L, M(2*j-1:2*j,2*l-1:2*l) = Bm(:,:,l-j+R+1); end, end
    M = (M-M.')/2;
end

function [Gd,mz] = tdl_corr(Bm,R)
    M = tdl_window(Bm,R+1); mz = real(majorana_string_ev([1 2],-1i,M));
    Gd = zeros(3,3,R+1); Gd(:,:,1) = diag([1 1 1-mz^2]);
    for r = 1:R
        for a = 1:2
            for b = 1:2
                Gd(a,b,r+1) = real(majorana_string_ev([a, 1:2*r, 2*r+b],(-1i)^r,M));
            end
        end
        Gd(3,3,r+1) = real(majorana_string_ev([1 2 2*r+1 2*r+2],-1,M)) - mz^2;
    end
end

function [Gd,mz,R] = tdl_corr_adaptive(mu,k)
    R = 40; Rmax = 256;
    while true
        [Gd,mz] = tdl_corr(tdl_blocks(mu,k,R),R);
        if max(max(abs(Gd(:,:,end)))) < 1e-7 || R >= Rmax, return; end
        R = min(2*R,Rmax);
    end
end

function [w,s,qa] = theorem1_from_corr(Gd)
    R = size(Gd,3)-1; q = linspace(0,pi,2001); S = zeros(3,numel(q));
    for a = 1:3
        c = squeeze(Gd(a,a,:)).';
        S(a,:) = c(1) + 2*c(2:end)*cos((1:R)'*q);
    end
    [s,ia] = min(S,[],2); qa = q(ia);
    w = min([sum(s)/2, s(1)+s(2), s(2)+s(3), s(3)+s(1)]);
end

function w = cmc_ti_sdp(Gd,nq)
    R = size(Gd,3)-1; q = 2*pi*(0:nq-1)/nq; Gq = zeros(3,3,nq);
    for n = 1:nq
        X = Gd(:,:,1);
        for r = 1:R, X = X + Gd(:,:,r+1)*exp(1i*q(n)*r) + Gd(:,:,r+1).'*exp(-1i*q(n)*r); end
        Gq(:,:,n) = (X+X')/2;
    end
    try
        cvx_begin quiet
            variable Z(3,3,nq) complex
            variable lam
            Z0 = zeros(3); obj = 0;
            for n = 1:nq
                Z(:,:,n) == hermitian_semidefinite(3);
                Z0 = Z0 + Z(:,:,n)/nq;
                obj = obj + real(trace(Z(:,:,n)*Gq(:,:,n)))/nq;
            end
            minimize(obj)
            subject to
                lam*eye(3) - Z0 == hermitian_semidefinite(3);
                real(trace(Z0)) - lam >= 1;
        cvx_end
        w = cvx_optval;
    catch err
        warning('cmc_ti_sdp failed (%s); using the structure-factor witness.',err.message);
        w = NaN;
    end
end

function m = pairwise_min_tdl(Gd,mz,rmax)
    s = {eye(2),[0 1;1 0],[0 -1i;1i 0],[1 0;0 -1]}; m = inf;
    for r = 1:min(rmax,size(Gd,3)-1)
        c = zeros(4); c(1,1) = 1; c(4,1) = mz; c(1,4) = mz;
        c(2:3,2:3) = Gd(1:2,1:2,r+1); c(4,4) = Gd(3,3,r+1) + mz^2;
        rho = zeros(4); for a = 0:3, for b = 0:3, rho = rho + c(a+1,b+1)*kron(s{a+1},s{b+1})/4; end, end
        T = reshape(permute(reshape(rho,[2 2 2 2]),[3 2 1 4]),4,4);
        m = min(m,min(real(eig((T+T')/2))));
    end
end

function rho = block_rdm(M,L)
    P = {eye(2),[0 1;1 0],[0 -1i;1i 0],[1 0;0 -1]}; rho = zeros(2^L);
    for code = 0:4^L-1
        lab = mod(floor(code./4.^(L-1:-1:0)),4); idx = []; ph = 1;
        for j = 1:L
            a = lab(j);
            if a == 3, idx = [idx 2*j-1 2*j]; ph = ph*(-1i);
            elseif a > 0
                for l = 1:j-1, idx = [idx 2*l-1 2*l]; ph = ph*(-1i); end
                idx = [idx 2*j-2+a];
            end
        end
        if isempty(idx), ev = 1; else, ev = majorana_string_ev(idx,ph,M); end
        if abs(ev) < 1e-14, continue; end
        Op = 1; for j = 1:L, Op = kron(Op,P{lab(j)+1}); end
        rho = rho + ev*Op;
    end
    rho = rho/2^L;
end

function m = block_ppt_min(M,L)
    rho = block_rdm(M,L); m = inf; cuts = all_cuts(L);
    for c = 1:numel(cuts), m = min(m,min(real(eig(pt_sites(rho,cuts{c},L))))); end
end

function tdl_selftest()
    J=1; gi=1; hi=0.3; gf=1; hf=2; Ti=1; t=1.5; L=8; N=120;
    S_i = tdl_setup(J,gi,hi,2048); S_f = tdl_setup(J,gf,hf,2048);
    Bm = tdl_blocks(tdl_mu_evolve(tdl_mu_gibbs(S_i.A,Ti),S_f.A,t),S_f.k,L);
    M = majorana_evolve(majorana_M_gibbs(majorana_A(N,J,gi,hi),Ti),majorana_A(N,J,gf,hf),t);
    s = N/2-L/2; I = 2*s+1:2*(s+L);
    d = max(max(abs(spin_covariance_ff(tdl_window(Bm,L),L) - spin_covariance_ff(M(I,I),L))));
    fprintf('Self-test infinite chain vs open chain N=%d (centre window): max|dGamma| = %.1e\n',N,d);
end

function A = majorana_A(N,J,g,h)
    A = zeros(2*N);
    for j = 1:N
        A(2*j-1,2*j) = A(2*j-1,2*j) + 2*h;
        if j < N
            A(2*j,2*j+1) = A(2*j,2*j+1) + J*(1+g);
            A(2*j-1,2*j+2) = A(2*j-1,2*j+2) - J*(1-g);
        end
    end
    A = A - A.';
end

function M = majorana_M_gibbs(A,T)
    X = 1i*A; X = (X+X')/2; [U,L] = eig(X); lam = real(diag(L));
    if T == 0, f = sign(lam); else, f = tanh(lam/(2*T)); end
    M = real(1i*(U*diag(f)*U')); M = (M-M.')/2;
end

function M = majorana_evolve(M0,A,t)
    O = expm(A*t); M = O*M0*O.'; M = (M-M.')/2;
end

function Mg = majorana_gge(M0,A)
    X = 1i*A; X = (X+X')/2; [U,L] = eig(X); lam = real(diag(L));
    D = abs(lam - lam.') < 1e-9;
    Mg = real(U*((U'*M0*U).*D)*U'); Mg = (Mg-Mg.')/2;
end

function [G,mz,C2] = spin_covariance_ff(M,N)
    mz = zeros(N,1);
    for j = 1:N, mz(j) = real(majorana_string_ev([2*j-1 2*j],-1i,M)); end
    G = zeros(3*N); C2 = zeros(N,N,3,3);
    for i = 1:N
        G(3*i-2:3*i,3*i-2:3*i) = diag([1 1 1-mz(i)^2]);
        for j = i+1:N
            mid = 2*i-1:2*j-2; ph = (-1i)^(j-i);
            for a = 1:2
                for b = 1:2
                    v = real(majorana_string_ev([2*i-2+a, mid, 2*j-2+b],ph,M));
                    C2(i,j,a,b) = v; G(3*i-3+a,3*j-3+b) = v; G(3*j-3+b,3*i-3+a) = v;
                end
            end
            v = real(majorana_string_ev([2*i-1 2*i 2*j-1 2*j],-1,M));
            C2(i,j,3,3) = v; G(3*i,3*j) = v - mz(i)*mz(j); G(3*j,3*i) = G(3*i,3*j);
        end
    end
end

function v = majorana_string_ev(idx,phase,M)
    sgn = 1; changed = true;
    while changed
        changed = false; p = 1;
        while p < numel(idx)
            if idx(p) > idx(p+1)
                idx([p p+1]) = idx([p+1 p]); sgn = -sgn; changed = true; p = p+1;
            elseif idx(p) == idx(p+1)
                idx(p:p+1) = []; changed = true; break;
            else
                p = p+1;
            end
        end
    end
    if isempty(idx), v = phase*sgn; return; end
    if mod(numel(idx),2) == 1, v = 0; return; end
    C = -1i*M(idx,idx); U = triu(C,1);
    v = phase*sgn*pfaffian_pr(U - U.');
end

function pf = pfaffian_pr(K)
    n = size(K,1); pf = 1;
    if mod(n,2) == 1, pf = 0; return; end
    for k = 1:2:n-1
        [~,p] = max(abs(K(k+1:n,k))); kp = k+p;
        if kp ~= k+1
            K([k+1 kp],:) = K([kp k+1],:); K(:,[k+1 kp]) = K(:,[kp k+1]); pf = -pf;
        end
        if K(k+1,k) == 0, pf = 0; return; end
        pf = pf*K(k,k+1);
        if k+2 <= n
            tau = K(k,k+2:n)/K(k,k+1); col = K(k+2:n,k+1);
            K(k+2:n,k+2:n) = K(k+2:n,k+2:n) + tau.'*col.' - col*tau;
        end
    end
end

function lmin = pairwise_min_ff(mz,C2,N)
    s = {eye(2),[0 1;1 0],[0 -1i;1i 0],[1 0;0 -1]}; lmin = inf;
    for i = 1:N-1
        for j = i+1:N
            c = zeros(4); c(1,1) = 1; c(4,1) = mz(i); c(1,4) = mz(j);
            c(2:4,2:4) = reshape(C2(i,j,:,:),3,3); r = zeros(4);
            for a = 0:3, for b = 0:3
                if c(a+1,b+1) ~= 0, r = r + c(a+1,b+1)*kron(s{a+1},s{b+1})/4; end
            end, end
            T = reshape(permute(reshape(r,[2 2 2 2]),[3 2 1 4]),4,4);
            lmin = min(lmin,min(real(eig((T+T')/2))));
        end
    end
end

function R = pt_sites(rho,sites,N)
    T = reshape(full(rho),2*ones(1,2*N)); p = 1:2*N;
    for s = sites, d = N-s+1; p([d d+N]) = p([d+N d]); end
    R = reshape(permute(T,p),2^N,2^N); R = (R+R')/2;
end

function cuts = all_cuts(N)
    cuts = {};
    for k = 1:floor(N/2)
        C = nchoosek(1:N,k);
        for r = 1:size(C,1)
            if 2*k == N && C(r,1) ~= 1, continue; end
            cuts{end+1} = C(r,:);
        end
    end
end

function x0 = crossing(x,y)
    x0 = NaN; k = find(y(1:end-1) < 1 & y(2:end) >= 1, 1, 'first');
    if ~isempty(k), x0 = x(k) + (1-y(k))*(x(k+1)-x(k))/(y(k+1)-y(k)); end
end

function H = local_H(N,J,g,h,op,P)
    H = sparse(2^N,2^N);
    for j = 1:N-1
        H = H - J*(1+g)/2*op(P{1},j)*op(P{1},j+1) - J*(1-g)/2*op(P{2},j)*op(P{2},j+1);
    end
    for j = 1:N, H = H - h*op(P{3},j); end
end

function ff_selftest()
    N = 6; J = 1; gi = 1.0; hi = 1.5; gf = 0.6; hf = 0.4; Ti = 0.8; t = 1.3;
    sx = sparse([0 1;1 0]); sy = sparse([0 -1i;1i 0]); sz = sparse([1 0;0 -1]); P = {sx,sy,sz};
    op = @(o,j) kron(kron(speye(2^(j-1)),o),speye(2^(N-j)));
    [V,D] = eig(full(local_H(N,J,gi,hi,op,P))); e = real(diag(D)); p = exp(-(e-min(e))/Ti); p = p/sum(p);
    [Vf,Df] = eig(full(local_H(N,J,gf,hf,op,P))); U = Vf*diag(exp(-1i*real(diag(Df))*t))*Vf';
    rho = U*(V*diag(p)*V')*U'; G_ed = zeros(3*N); m = zeros(N,3);
    for j = 1:N, for a = 1:3, m(j,a) = real(trace(rho*op(P{a},j))); end, end
    for j = 1:N, for a = 1:3, for k = 1:N, for b = 1:3
        X = op(P{a},j); Y = op(P{b},k);
        G_ed(3*(j-1)+a,3*(k-1)+b) = 0.5*real(trace(rho*(X*Y+Y*X))) - m(j,a)*m(k,b);
    end, end, end, end
    G_ff = spin_covariance_ff(majorana_evolve(majorana_M_gibbs(majorana_A(N,J,gi,hi),Ti),majorana_A(N,J,gf,hf),t),N);
    fprintf('Self-test open chain N=6 vs exact diagonalization: max|dGamma| = %.1e\n',max(abs(G_ed(:)-G_ff(:))));
end

function shade(ax,x,y,c)
    fill(ax,[x(1) x(2) x(2) x(1)],[y(1) y(1) y(2) y(2)],c,'EdgeColor','none','HandleVisibility','off');
end

function fmt(ax,FONT,FT,xl,yl)
    set(ax,'FontName',FONT,'FontSize',FT,'TickDir','out','Layer','top','XLim',xl,'YLim',yl); grid(ax,'on');
end

function save_fig(fig,outdir,name)
    exportgraphics(fig,fullfile(outdir,[name '.pdf']),'ContentType','vector');
    exportgraphics(fig,fullfile(outdir,[name '.png']),'Resolution',300);
end
