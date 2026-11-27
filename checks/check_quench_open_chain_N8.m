% CHECK_QUENCH_OPEN_CHAIN_N8  Open chain of 8 spins after the quench: all 127 cuts and the CMC (App. G).
% Part 2 needs CVX. Several minutes. Expected: NPT at all 40 times t = 0.5, 1, ..., 20;
% CMC undetected only at t = 0, 4.5 and 5.

%% Setup
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));
has_cvx = exist('cvx_begin', 'file') > 0;

N = 8;  gamma = 1.0;  h_i = 0.3;  h_f = 2.0;  T_i = 1.0;
ops = local_pauli_operators(N);
rho0 = gibbs_state(xy_hamiltonian(N, gamma, h_i, false), T_i);
[Vf, Ef] = eig(full(xy_hamiltonian(N, gamma, h_f, false)));
Ef = diag(Ef);
rho_at = @(t) (Vf * diag(exp(-1i * Ef * t)) * Vf') * rho0 * (Vf * diag(exp(-1i * Ef * t)) * Vf')';

cuts = {};
for m = 1:N-1
    subsets = nchoosek(2:N, m);
    for r = 1:size(subsets, 1)
        cuts{end+1} = subsets(r, :);
    end
end
fprintf('%d bipartitions\n', numel(cuts));

%% Part 1: the initial state is PPT across every bipartition
values = zeros(1, numel(cuts));
for n = 1:numel(cuts)
    values(n) = min(real(eig(partial_transpose(rho0, cuts{n}, N))));
end
fprintf('t = 0: smallest partial-transpose eigenvalue over all cuts = %.3e; cuts with a negative eigenvalue: %d\n', ...
    min(values), sum(values < -1e-12));
fprintf('       smallest eigenvalue of rho(0) itself = %.3e\n', min(real(eig(rho0))));

%% Part 2: partial transposes and CMC after the quench
grids = {0:0.5:20, 0:0.25:10};
for g = 1:numel(grids)
    times = grids{g};
    pt_min = zeros(size(times));
    omega = nan(size(times));
    for n = 1:numel(times)
        rho = rho_at(times(n));
        rho = (rho + rho') / 2;
        m = Inf;
        for c = 1:numel(cuts)
            m = min(m, min(real(eig(partial_transpose(rho, cuts{c}, N)))));
        end
        pt_min(n) = m;
        if has_cvx
            omega(n) = cmc_sdp(covariance_matrix(rho, ops));
        end
    end
    fprintf('\nt = %g .. %g in steps of %g\n', times(1), times(end), times(2) - times(1));
    fprintf('  t = 0: smallest PT eigenvalue over all cuts = %.4f, omega_CMC = %.4f\n', pt_min(1), omega(1));
    fprintf('  t > 0: PT negative at %d of %d times\n', sum(pt_min(2:end) < 0), numel(times) - 1);
    if has_cvx
        fprintf('  CMC detects (omega < 1) at %d of %d times; undetected times: %s\n', ...
            sum(omega < 1), numel(times), mat2str(times(omega >= 1)));
        fprintf('  omega range for t > 0: %.3f .. %.3f\n', min(omega(2:end)), max(omega(2:end)));
    else
        fprintf('  CMC skipped (CVX not found)\n');
    end
end
