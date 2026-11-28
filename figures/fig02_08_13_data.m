% FIG02_08_13_DATA  Moment hierarchy (k = 2), CMC, PPT witness and PPT mixer versus T, N = 4 (Figs. 2, 8, 13).
% Needs CVX (MOSEK recommended). The PPT mixer needs fpptwit and YALMIP on the path; it is skipped otherwise.
% only_missing = true computes only the (gamma, h) points that are not yet in data/N4_criteria_vs_T_<BC>.mat
% and appends them; only_missing = false recomputes everything.

repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));

only_missing = true;
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
nT = numel(T_list);

for BC = {'open', 'periodic'}
    file = fullfile(repo, 'data', sprintf('N%d_criteria_vs_T_%s.mat', N, BC{1}));
    todo = 1:size(cases, 1);
    old = [];
    if only_missing && exist(file, 'file')
        old = load(file, 'settings', 'results');
        if numel(old.settings.T_list) ~= nT || max(abs(old.settings.T_list(:) - T_list(:))) > 1e-12
            error('The temperature grid of %s differs; set only_missing = false.', file);
        end
        have = false(1, size(cases, 1));
        for ic = 1:size(cases, 1)
            have(ic) = any(all(abs(old.settings.parameterCases - cases(ic, :)) < 1e-9, 2));
        end
        todo = find(~have);
    end
    if isempty(todo)
        fprintf('%s: all points already in %s\n', BC{1}, file);
        continue
    end
    ENTSDP = zeros(numel(todo), nT);
    ENTCMC = zeros(numel(todo), nT);
    ENTPPT = zeros(numel(todo), nT);
    FPPTwitStrength = zeros(numel(todo), nT);
    FPPTwitRaw = nan(numel(todo), nT);
    for n = 1:numel(todo)
        gh = cases(todo(n), :);
        stopped = false(1, 4);
        for it = 1:nT
            rho = gibbs_state(xy_hamiltonian(N, gh(1), gh(2), strcmp(BC{1}, 'periodic')), T_list(it));
            if ~stopped(1)
                ENTSDP(n, it) = max(0, -moment_hierarchy_value(rho, N, k, solver));
                stopped(1) = ENTSDP(n, it) <= detectTol;
            end
            if ~stopped(2)
                ENTCMC(n, it) = max(0, 1 - cmc_sdp(covariance_matrix(rho, ops)));
                stopped(2) = ENTCMC(n, it) <= detectTol;
            end
            if ~stopped(3)
                ENTPPT(n, it) = max(0, -ppt_negative_part(rho, cut_A, N));
                stopped(3) = ENTPPT(n, it) <= detectTol;
            end
            if hasFpptwit && ~stopped(4)
                FPPTwitRaw(n, it) = real(fpptwit(rho));
                FPPTwitStrength(n, it) = max(0, -FPPTwitRaw(n, it));
                stopped(4) = FPPTwitStrength(n, it) <= detectTol;
            end
        end
        fprintf('%s, (gamma, h) = (%g, %g) done\n', BC{1}, gh(1), gh(2));
    end
    newCases = cases(todo, :);
    if ~isempty(old)
        R0 = old.results(1);
        newCases = [old.settings.parameterCases; newCases];
        ENTSDP = [R0.ENTSDP; ENTSDP];
        ENTCMC = [R0.ENTCMC; ENTCMC];
        ENTPPT = [R0.ENTPPT; ENTPPT];
        FPPTwitStrength = [R0.FPPTwitStrength; FPPTwitStrength];
        FPPTwitRaw = [R0.FPPTwitRaw; FPPTwitRaw];
    end
    settings = struct('Nq', N, 'BC', BC{1}, 'kDual', k, 'solverName', solver, ...
        'parameterCases', newCases, 'T_list', T_list, 'detectTol', detectTol, 'B_sites', cut_A);
    results = struct('ENTSDP', ENTSDP, 'ENTCMC', ENTCMC, 'ENTPPT', ENTPPT, ...
        'FPPTwitStrength', FPPTwitStrength, 'FPPTwitRaw', FPPTwitRaw);
    save(file, 'settings', 'results', '-v7.3');
    fprintf('%s: saved %s\n', BC{1}, file);
end
