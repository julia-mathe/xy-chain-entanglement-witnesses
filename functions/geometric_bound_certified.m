function [bound, lambda_best, eta_best] = geometric_bound_certified(omega, N, C, ells)
% GEOMETRIC_BOUND_CERTIFIED  Certified lower bound on E_G from Eq. (36) with the PPT relaxation (App. H).

require_cvx();
w_measured = omega - 1;
w = 2*ells(:).*(ells(:) + 1)/N - 1;
number_of_conditions = size(C, 1);

cvx_begin quiet
    variable lambda
    variable eta
    variable u(number_of_conditions) nonnegative
    maximize( 1 + lambda * w_measured - eta )
    subject to
        lambda <= 0;
        inv_pos(eta - lambda * w) + C' * u <= 1;
cvx_end

bound = cvx_optval;
lambda_best = lambda;
eta_best = eta;
end
