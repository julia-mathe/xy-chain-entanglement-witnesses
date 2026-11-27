% CHECK_FULLY_PPT_OPEN_CHAIN_N4  CMC versus PPT (all cuts) in the (gamma, h) plane, N = 4, OBC, T = 1.47
% (Sec. 3, App. D).
% Needs CVX. 10-30 min. Expected: detected 701, NPT 760, detected and fully PPT 17
% (|gamma| <= 0.1, |h| <= 0.3, omega_CMC between 0.994 and 0.998).

%% Setup
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));

N = 4;
T = 1.47;
ops = local_pauli_operators(N);
cuts = {1, 2, 3, 4, [1 2], [1 3], [1 4]};
gammas = linspace(-1.5, 1.5, 31);
fields = linspace(-2.5, 2.5, 51);

%% Scan
results = zeros(numel(gammas) * numel(fields), 5);
row = 0;
tic;
for gamma = gammas
    for h = fields
        rho = gibbs_state(xy_hamiltonian(N, gamma, h, false), T);
        Gamma = covariance_matrix(rho, ops);
        omega = cmc_sdp(Gamma);
        smallest = zeros(1, numel(cuts));
        for c = 1:numel(cuts)
            smallest(c) = min(eig(partial_transpose(rho, cuts{c}, N)));
        end
        [worst, which_cut] = min(smallest);
        row = row + 1;
        results(row, :) = [gamma, h, omega, worst, which_cut];
    end
    fprintf('gamma = %+.1f done (%.0f s)\n', gamma, toc);
end

%% Summary
detected  = results(:, 3) < 1 - 1e-6;
fully_ppt = results(:, 4) > -1e-10;
fprintf('\ndetected by the CMC: %d,  NPT: %d,  detected and fully PPT: %d,  NPT and not detected: %d\n', ...
        sum(detected), sum(~fully_ppt), sum(detected & fully_ppt), sum(~fully_ppt & ~detected));

fprintf('\nPoints detected by the CMC although PPT across every cut:\n');
fprintf('%8s %8s %10s %14s   %s\n', 'gamma', 'h', 'omega_CMC', 'min PT eig.', 'transposed spins');
for r = find(detected & fully_ppt)'
    fprintf('%8.2f %8.2f %10.5f %14.2e   %s\n', results(r, 1), results(r, 2), results(r, 3), ...
            results(r, 4), mat2str(cuts{results(r, 5)}));
end

if ~exist(fullfile(repo, 'output'), 'dir'), mkdir(fullfile(repo, 'output')); end
save(fullfile(repo, 'output', 'fully_ppt_scan_N4_OBC.mat'), 'results', 'gammas', 'fields', 'T', 'N', '-v7');
