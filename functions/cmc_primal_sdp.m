function [s, X] = cmc_primal_sdp(Gamma)
% CMC_PRIMAL_SDP  Largest s with Gamma - (+)_j (1 - X_j) - s*1 >= 0, Tr X_j = 1 (Eq. (9)).

require_cvx();
N = size(Gamma, 1) / 3;

cvx_begin quiet
    variable s
    variable Xall(3,3*N)
    expression LocalNoise(3*N,3*N)
    maximize( s )
    subject to
        for j = 1:N
            block = 3*j-2 : 3*j;
            Xall(:, block) == semidefinite(3);
            trace(Xall(:, block)) == 1;
            LocalNoise(block, block) = eye(3) - Xall(:, block);
        end
        Gamma - LocalNoise - s * eye(3*N) == semidefinite(3*N);
cvx_end

X = zeros(3, 3, N);
for j = 1:N
    X(:, :, j) = Xall(:, 3*j-2 : 3*j);
end
end
