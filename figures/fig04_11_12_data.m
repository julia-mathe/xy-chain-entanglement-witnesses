% FIG04_11_12_DATA  Gamma and Z* for N = 10, T = 0.1 (data of Figs. 4, 11 and 12).
% Needs CVX. About a minute. Writes data/matrix_figures_N10_T0.1.mat.

%% Setup
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));
require_cvx();
if ~exist(fullfile(repo, 'data'), 'dir'), mkdir(fullfile(repo, 'data')); end

N = 10;  T = 0.1;
ops = local_pauli_operators(N);
definitions = { 'PBC_Ising', 1.0, 1.0, true
                'PBC_XX',    0.0, 0.0, true
                'OBC_Ising', 1.0, 1.0, false
                'OBC_XX',    0.0, 0.0, false
                'PBC_nc',    0.5, 0.5, true
                'OBC_nc',    0.5, 0.5, false };

cases = struct('name', {}, 'gamma', {}, 'h', {}, 'periodic', {}, 'Gamma', {}, 'Z', {}, 'omega', {}, 'b', {});
for n = 1:size(definitions, 1)
    tic;
    [name, gamma, h, periodic] = definitions{n, :};
    rho = gibbs_state(xy_hamiltonian(N, gamma, h, periodic), T);
    Gamma = covariance_matrix(rho, ops);
    [omega, Z] = cmc_sdp(Gamma);
    Z = (Z + Z') / 2;
    cases(n) = struct('name', name, 'gamma', gamma, 'h', h, 'periodic', periodic, ...
                      'Gamma', Gamma, 'Z', Z, 'omega', omega, 'b', witness_bound_b(Z));
    ev = sort(eig(Z), 'descend');
    fprintf('%-10s omega_SDP = %.5f  b = %.6f  largest eig(Z*) = [%s]', name, omega, cases(n).b, sprintf(' %.4f', ev(1:4)));
    if periodic
        [omega_th, s, qmin, best] = cmc_theorem1_ring(Gamma);
        fprintf('  Theorem 1 = %.5f (%s), s = [%.3f %.3f %.3f], q/pi = [%.3f %.3f %.3f]', omega_th, best, s, qmin/pi);
    end
    fprintf('  [%.0f s]\n', toc);
end

N_big = 200;  L = 10;
G = majorana_G_open(N_big, 1.0, 1.0, T);
Gamma = covariance_window_open(G, N_big/2 - L/2 + 1, L);
[omega, Z] = cmc_sdp(Gamma);
Z = (Z + Z') / 2;
cases(7) = struct('name', 'window', 'gamma', 1.0, 'h', 1.0, 'periodic', false, ...
                  'Gamma', Gamma, 'Z', Z, 'omega', omega, 'b', witness_bound_b(Z));
fprintf('%-10s omega_SDP = %.5f  b = %.6f\n', 'window', omega, cases(7).b);

save(fullfile(repo, 'data', 'matrix_figures_N10_T0.1.mat'), 'cases', 'N', 'T', '-v7');
fprintf('Saved data/matrix_figures_N10_T0.1.mat\n');
