% CHECK_FULLY_PPT_EXAMPLE_N4  Fully PPT thermal state detected by the CMC: XX point, OBC, N = 4, T = 1.47
% (Sec. 3, App. D).
% Needs CVX. Expected: all seven partial transposes positive (smallest 5.1e-4), omega_CMC = 0.99437,
% primal value s < 0.

%% Setup
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));

N = 4;
T = 1.47;
ops = local_pauli_operators(N);
rho = gibbs_state(xy_hamiltonian(N, 0, 0, false), T);
cuts = {1, 2, 3, 4, [1 2], [1 3], [1 4]};

%% 1. Partial transposes (two independent implementations)
fprintf('Smallest eigenvalue of the partial transpose:\n');
for c = 1:numel(cuts)
    fast = min(eig(partial_transpose(rho, cuts{c}, N)));
    loops = min(eig(partial_transpose_by_loops(rho, cuts{c}, N)));
    fprintf('  transposed spins %-6s:  %.6e   (loop version: %.6e)\n', mat2str(cuts{c}), fast, loops);
end

%% 2. Optimal witness from the SDP
Gamma = covariance_matrix(rho, ops);
[omega, Z] = cmc_sdp(Gamma);
Z = (Z + Z') / 2;
fprintf('\nomega_CMC (SDP) = %.6f\n', omega);
fprintf('smallest eigenvalue of Z = %.2e\n', min(eig(Z)));
fprintf('Tr(Z*Gamma) = %.8f,  b(Z) = %.8f,  ratio = %.8f\n', trace(Z*Gamma), witness_bound_b(Z), ...
        trace(Z*Gamma) / witness_bound_b(Z));
[U, e] = eig(Z);
Z_clean = U * diag(max(diag(e), 0)) * U';
fprintf('cleaned Z: ratio Tr(Z*Gamma)/b(Z) = %.8f\n', trace(Z_clean*Gamma) / witness_bound_b(Z_clean));
eigenvalues = sort(eig(Z), 'descend');
fprintf('largest eigenvalues of Z: %s\n', sprintf('%.4f  ', eigenvalues(1:6)));

%% 3. Primal form of the criterion
s = cmc_primal_sdp(Gamma);
fprintf('\nprimal problem: largest s = %.6f  (negative => covariance matrix not of separable form)\n', s);

%% 4. Purity
fprintf('purity Tr(rho^2) = %.8f\n', real(trace(rho * rho)));
