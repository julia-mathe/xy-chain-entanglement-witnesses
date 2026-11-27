function Spi = staggered_structure_factor_euler(C, depth)
% STAGGERED_STRUCTURE_FACTOR_EULER  S(pi) by repeated averaging of partial sums (App. A).

if nargin < 2
    depth = 6;
end

r = 0:numel(C) - 1;
terms = 2 * (-1).^r .* C;
terms(1) = C(1);
P = cumsum(terms);

for level = 1:depth
    P = (P(2:end) + P(1:end-1)) / 2;
end
Spi = P(end);
end
