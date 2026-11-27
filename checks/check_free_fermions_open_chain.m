% CHECK_FREE_FERMIONS_OPEN_CHAIN  Free fermions versus exact diagonalization (App. A);
% window of 10 spins in an open chain of 200 (Sec. 4.1, Fig. 11).
% Part 2 needs CVX. Expected: differences ~1e-15; window omega_CMC = 0.59113.

%% Setup
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));

%% Part 1: exact diagonalization versus free fermions, open chain of 8 spins
N = 8;
ops = local_pauli_operators(N);
parameters = [1.0, 1.0, 0.5
              0.3, 0.7, 0.4];
for n = 1:size(parameters, 1)
    gamma = parameters(n, 1);   h = parameters(n, 2);   T = parameters(n, 3);
    rho = gibbs_state(xy_hamiltonian(N, gamma, h, false), T);
    Gamma_ed = covariance_matrix(rho, ops);
    G = majorana_G_open(N, gamma, h, T);
    Gamma_ff = covariance_window_open(G, 1, N);
    fprintf('gamma = %.1f, h = %.1f, T = %.1f:  largest |Gamma_ED - Gamma_FF| = %.1e\n', ...
            gamma, h, T, max(abs(Gamma_ed(:) - Gamma_ff(:))));
end

%% Part 2: 10 central spins of an open chain of 200 spins
N_big = 200;   L = 10;
G = majorana_G_open(N_big, 1.0, 1.0, 0.1);
first_site = (N_big - L)/2 + 1;
Gamma_window = covariance_window_open(G, first_site, L);
omega_window = cmc_sdp(Gamma_window);
fprintf('\nWindow of %d spins in an open chain of %d (Ising point, T = 0.1): omega_CMC = %.5f\n', ...
        L, N_big, omega_window);
