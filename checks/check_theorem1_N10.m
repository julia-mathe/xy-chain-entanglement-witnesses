% CHECK_THEOREM1_N10  SDP versus Theorem 1, N = 10, T = 0.1 (Sec. 4.2, Tables 1 and 2, Fig. 4).
% Needs CVX. About a minute.
% Expected: Ising PBC 0.58269 (= Theorem 1), OBC 0.45972 (three-mode 0.47446);
%           XX PBC 0.15371 (= Theorem 1), OBC 0.13980 (three-mode 0.14375);
%           (0.5,0.5) PBC 0.77770 (q_x = 0.8 pi, Z* eigenvalues 1, 1/2, 1/2; three-mode 0.78759),
%           OBC 0.77218 (three-mode 0.79437).

%% Setup: make the building blocks in the folder "functions" available
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));

%% Parameters
N = 10;
T = 0.1;
ops = local_pauli_operators(N);

points = { 1,   1,   'Ising point (1, 1)'
           0,   0,   'XX point (0, 0)'
           0.5, 0.5, 'noncritical point (0.5, 0.5)' };

%% Main loop over the points and the boundary conditions
for n = 1:size(points, 1)
    gamma = points{n, 1};
    h     = points{n, 2};
    name  = points{n, 3};

    for periodic = [true, false]
        H = xy_hamiltonian(N, gamma, h, periodic);
        rho = gibbs_state(H, T);
        [Gamma, m] = covariance_matrix(rho, ops);
        [omega_sdp, Z] = cmc_sdp(Gamma);

        if periodic
            boundary = 'PBC';
        else
            boundary = 'OBC';
        end
        fprintf('\n%s, %s, N = %d, T = %.2f\n', name, boundary, N, T);
        fprintf('  omega_CMC from the SDP             = %.5f\n', omega_sdp);

        if periodic
            [omega_th, s, qmin, best] = cmc_theorem1_ring(Gamma);
            [~, ~, candidates] = theorem1_omega(s);
            fprintf('  omega_CMC from Theorem 1           = %.5f   (smallest candidate: %s)\n', omega_th, best);
            fprintf('  candidates [half-sum, xy, yz, xz]  = [%.5f %.5f %.5f %.5f]\n', candidates);
            fprintf('  minima [s_x s_y s_z]               = [%.4f %.4f %.4f]\n', s);
            fprintf('  at momenta q/pi                    = [%.3f %.3f %.3f]\n', qmin / pi);
        end

        fprintf('  XX-type witness    (q = pi, pi, 0) = %.5f\n', three_mode_witness(Gamma, [pi, pi, 0]));
        fprintf('  Ising-type witness (q = pi, 0, pi) = %.5f\n', three_mode_witness(Gamma, [pi, 0, pi]));
        fprintf('  mean <sz> = %.5f,  max |<sx>| = %.1e,  max |<sy>| = %.1e\n', ...
                mean(m(3:3:end)), max(abs(m(1:3:end))), max(abs(m(2:3:end))));
        eigenvalues = sort(eig((Z + Z')/2), 'descend');
        fprintf('  largest eigenvalues of Z*          = %s\n', sprintf('%.4f  ', eigenvalues(1:6)));
    end
end
