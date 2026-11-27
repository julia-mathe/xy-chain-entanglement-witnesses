% CHECK_MOMENT_HIERARCHY_N4  First level of the moment hierarchy with one- and two-spin data equals the CMC
% (App. C).
% Needs CVX. Expected: equal noise thresholds, e.g. 0.86611 for both at the XX point, T = 1.47, PBC.

%% Setup
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));

N = 4;
k = 2;
ops = local_pauli_operators(N);

cases = { 0, 0, 1.47, true
          1, 1, 1.47, true
          0, 0, 1.00, true
          1, 1, 1.00, true
          0, 0, 1.47, false };

fprintf('%-26s %10s %10s\n', '(gamma, h, T, boundary)', 'p_moment', 'p_CMC');
for n = 1:size(cases, 1)
    [gamma, h, T, periodic] = cases{n, :};
    rho = gibbs_state(xy_hamiltonian(N, gamma, h, periodic), T);
    [~, m, C] = covariance_matrix(rho, ops);

    p_moment = moment_hierarchy_visibility(m, C, N, k);
    p_cmc = cmc_visibility(m, C, N, periodic);

    if periodic
        boundary = 'periodic';
    else
        boundary = 'open';
    end
    fprintf('(%g, %g, %.2f, %-8s)     %10.5f %10.5f\n', gamma, h, T, boundary, p_moment, p_cmc);
end
