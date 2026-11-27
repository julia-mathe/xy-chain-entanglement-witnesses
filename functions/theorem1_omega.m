function [omega, best, candidates] = theorem1_omega(s)
% THEOREM1_OMEGA  Right-hand side of Eq. (25) from the minimal structure factors s.

candidates = [ (s(1) + s(2) + s(3))/2,   s(1) + s(2),   s(2) + s(3),   s(1) + s(3) ];
names      = { 'half-sum',               'xy',          'yz',          'xz'        };

[omega, k] = min(candidates);
best = names{k};
end
