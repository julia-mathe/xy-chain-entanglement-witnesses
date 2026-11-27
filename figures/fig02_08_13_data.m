% FIG02_08_13_DATA  Moment hierarchy (k = 2), CMC, PPT witness and PPT mixer versus T, N = 4 (Figs. 2, 8, 13).
% Needs CVX (MOSEK recommended). The PPT mixer needs fpptwit and YALMIP on the path; it is skipped otherwise.
% Writes data/N4_criteria_vs_T_open.mat and data/N4_criteria_vs_T_periodic.mat.

repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));

N = 4;
k = 2;
solver = 'mosek';
cases = [0 0; 1 1; 0.6 0.2; 0.9 0.9; 1.3 1.2];
T_list = linspace(0.10, 2.60, 260);
cut_A = 2:2:N;
detectTol = 1e-6;
hasFpptwit = exist('fpptwit', 'file') == 2;
if ~hasFpptwit
    fprintf('fpptwit not found: PPT mixer (Fig. 13) skipped.\n');
end
ops = local_pauli_operators(N);

for BC = ["open", "periodic"]
    nC = size(cases, 1);
    nT = numel(T_list);
    ENTSDP = nan(nC, nT);
    ENTCMC = nan(nC, nT);
    ENTPPT = nan(nC, nT);
    FPPTwitStrength = nan(nC, nT);
    FPPTwitRaw = nan(nC, nT);
    for ic = 1:nC
        stopped = false(1, 4);
        for it = 1:nT
            rho = gibbs_state(xy_hamiltonian(N, cases(ic, 1), cases(ic, 2), BC == "periodic"), T_list(it));
            ENTSDP(ic, it) = 0;
            ENTCMC(ic, it) = 0;
            ENTPPT(ic, it) = 0;
            FPPTwitStrength(ic, it) = 0;
            if ~stopped(1)
                ENTSDP(ic, it) = max(0, -moment_hierarchy_value(rho, N, k, solver));
                stopped(1) = ENTSDP(ic, it) <= detectTol;
            end
            if ~stopped(2)
                ENTCMC(ic, it) = max(0, 1 - cmc_sdp(covariance_matrix(rho, ops)));
                stopped(2) = ENTCMC(ic, it) <= detectTol;
            end
            if ~stopped(3)
                ENTPPT(ic, it) = max(0, -ppt_negative_part(rho, cut_A, N));
                stopped(3) = ENTPPT(ic, it) <= detectTol;
            end
            if hasFpptwit && ~stopped(4)
                FPPTwitRaw(ic, it) = real(fpptwit(rho));
                FPPTwitStrength(ic, it) = max(0, -FPPTwitRaw(ic, it));
                stopped(4) = FPPTwitStrength(ic, it) <= detectTol;
            end
        end
        fprintf('%s, (gamma, h) = (%g, %g) done\n', BC, cases(ic, 1), cases(ic, 2));
    end
    settings = struct('Nq', N, 'BC', char(BC), 'kDual', k, 'solverName', solver, ...
        'parameterCases', cases, 'T_list', T_list, 'detectTol', detectTol, 'B_sites', cut_A);
    results = struct('ENTSDP', ENTSDP, 'ENTCMC', ENTCMC, 'ENTPPT', ENTPPT, ...
        'FPPTwitStrength', FPPTwitStrength, 'FPPTwitRaw', FPPTwitRaw);
    save(fullfile(repo, 'data', sprintf('N%d_criteria_vs_T_%s.mat', N, BC)), 'settings', 'results', '-v7.3');
end
