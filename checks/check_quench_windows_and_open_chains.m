% CHECK_QUENCH_WINDOWS_AND_OPEN_CHAINS  Steady state of windows of L spins and of open chains (Fig. 14(c)).
% Needs CVX. A few minutes. Expected: L = 2, 3, 4, 8, 12 give 1.0776, 0.9902, 0.9555, 0.9172, 0.9090;
% open chains N = 8, 12, 16, 20 give 0.9066, 0.9020, 0.9007, 0.9003.

%% Setup
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));
require_cvx();

gamma = 1.0;  h_i = 0.3;  h_f = 2.0;  T_i = 1.0;
Mg = majorana_state_symbol_tdl('steady', majorana_state_symbol_tdl('thermal', gamma, h_i, T_i), gamma, h_f);

%% Part 1: windows of the infinite-chain steady state
fprintf('windows of the infinite-chain GGE:\n');
for L = 2:12
    Mw = majorana_window_tdl(Mg, L);
    omega = cmc_sdp(covariance_from_majorana(Mw, 1:L));
    fprintf('  L = %2d: omega_CMC = %.4f\n', L, omega);
end

%% Part 2: the window L = 3 from its reduced density matrix
Mw = majorana_window_tdl(Mg, 3);
rho3 = reduced_density_matrix_gaussian(Mw, 1:3);
Gamma_rho = covariance_matrix(rho3, local_pauli_operators(3));
Gamma_ff = covariance_from_majorana(Mw, 1:3);
fprintf('\nrho_3: trace %.6f, smallest eigenvalue %.4f; max |Gamma(rho_3) - Gamma(Wick)| = %.1e; omega_CMC = %.5f\n', ...
    real(trace(rho3)), min(real(eig(rho3))), max(abs(Gamma_rho(:) - Gamma_ff(:))), cmc_sdp(Gamma_rho));

%% Part 3: open chains
fprintf('\nopen chains (steady state of the finite chain):\n');
for N = 8:2:20
    tic;
    M0 = majorana_thermal_state(majorana_hamiltonian_matrix(N, gamma, h_i), T_i);
    Mss = majorana_steady_state_open(M0, majorana_hamiltonian_matrix(N, gamma, h_f));
    omega = cmc_sdp(covariance_from_majorana(Mss, 1:N));
    fprintf('  N = %2d: steady-state omega_CMC = %.4f  [%.0f s]\n', N, omega, toc);
end
