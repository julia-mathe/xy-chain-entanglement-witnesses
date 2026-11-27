function [bound, c_best] = geometric_bound_large_N(omega)
% GEOMETRIC_BOUND_LARGE_N  Large-N limit of the Neel-type bound on E_G (App. H).

minus_bound = @(logc) -(1 - (omega + exp(logc)) .* exp(exp(logc)) .* expint(exp(logc)));
[logc_best, value] = fminbnd(minus_bound, -12, 6, optimset('TolX', 1e-12));
bound = -value;
c_best = exp(logc_best);
end
