function omega = cmc_block_sdp_tdl(c, nq)
% CMC_BLOCK_SDP_TDL  Optimal CMC value of a translation-invariant infinite chain, Eq. (58).

require_cvx();
if nargin < 2 || isempty(nq)
    nq = 64;
end
q = 2*pi*(0:nq-1) / nq;
G = covariance_fourier_tdl(c, q);

cvx_begin quiet
    variable Z(3, 3, nq) complex
    variable t
    Z0 = sum(Z, 3) / nq;
    objective = 0;
    for n = 1:nq
        objective = objective + real(trace(Z(:, :, n) * G(:, :, n)));
    end
    minimize( objective / nq )
    subject to
        for n = 1:nq
            Z(:, :, n) == hermitian_semidefinite(3);
        end
        t * eye(3) - Z0 == hermitian_semidefinite(3);
        real(trace(Z0)) - t >= 1;
cvx_end

if ~strcmp(cvx_status, 'Solved')
    warning('cmc_block_sdp_tdl: the solver returned the status "%s".', cvx_status);
end
omega = cvx_optval;
end
