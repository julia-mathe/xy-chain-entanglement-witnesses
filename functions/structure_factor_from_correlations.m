function S = structure_factor_from_correlations(C, q)
% STRUCTURE_FACTOR_FROM_CORRELATIONS  Structure factor of the infinite chain, Eq. (29).

r = 0:numel(C) - 1;
S = zeros(size(q));
for n = 1:numel(q)
    S(n) = C(1) + 2 * sum( C(2:end) .* cos(q(n) * r(2:end)) );
end
end
