% CHECK_THEOREM1_RANDOM_PARAMETERS  SDP versus Theorem 1 at 12 random (gamma, h, T), N = 6 (Sec. 4.2).
% Needs CVX. Expected: largest difference ~1e-8.

%% Setup
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));

N = 6;
number_of_points = 12;
ops = local_pauli_operators(N);
rng(1);

x = 1:3:3*N;   y = 2:3:3*N;   z = 3:3:3*N;

%% Loop over random points
largest_difference = 0;
fprintf('%8s %8s %8s | %10s %10s %10s | %9s  %s\n', 'gamma', 'h', 'T', 'SDP', 'Theorem 1', ...
        'offblocks', 'best', 'q/pi');
for trial = 1:number_of_points
    gamma = -1.5 + 3 * rand();
    h     = -2 + 4 * rand();
    T     = 10^(-1.5 + 2 * rand());

    rho = gibbs_state(xy_hamiltonian(N, gamma, h, true), T);
    Gamma = covariance_matrix(rho, ops);
    offblocks = max([max(max(abs(Gamma(x, y)))), max(max(abs(Gamma(x, z)))), max(max(abs(Gamma(y, z))))]);

    omega_sdp = cmc_sdp(Gamma);
    [omega_th, s, qmin, best] = cmc_theorem1_ring(Gamma);
    largest_difference = max(largest_difference, abs(omega_sdp - omega_th));

    fprintf('%+8.3f %+8.3f %8.3f | %10.5f %10.5f %10.1e | %9s  [%s]\n', gamma, h, T, omega_sdp, ...
            omega_th, offblocks, best, sprintf('%.3f ', qmin / pi));
end
fprintf('\nLargest difference between SDP and Theorem 1: %.1e\n', largest_difference);
