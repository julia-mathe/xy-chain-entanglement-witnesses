function [omega, Z] = cmc_sdp(Gamma)
% CMC_SDP  Optimal CMC value omega_CMC and witness Z, Eqs. (10)-(11).

require_cvx();
N = size(Gamma, 1) / 3;

cvx_begin sdp quiet
    variable Z(3*N,3*N) symmetric
    variable mu(N)
    minimize( trace(Z * Gamma) )
    subject to
        Z >= 0;
        bZ = 0;
        for j = 1:N
            block = 3*j-2 : 3*j;
            mu(j) * eye(3) >= Z(block, block);
            bZ = bZ + trace(Z(block, block)) - mu(j);
        end
        bZ >= 1;
cvx_end

if ~strcmp(cvx_status, 'Solved')
    warning('cmc_sdp: the solver returned the status "%s".', cvx_status);
end
omega = cvx_optval;
end
