function [ells, P, D, ell_of_state, V] = total_spin_projectors(N)
% TOTAL_SPIN_PROJECTORS  Projectors onto the total-spin sectors of N spins 1/2.

ops = local_pauli_operators(N);
Sx = 0;  Sy = 0;  Sz = 0;
for j = 1:N
    Sx = Sx + ops{j,1} / 2;
    Sy = Sy + ops{j,2} / 2;
    Sz = Sz + ops{j,3} / 2;
end
J2 = real(full(Sx*Sx + Sy*Sy + Sz*Sz));

[V, ev] = eig(J2);
ev = diag(ev);
ell_of_state = round(2 * (-1 + sqrt(1 + 4*ev)) / 2) / 2;

ells = unique(ell_of_state)';
P = cell(1, numel(ells));
D = zeros(1, numel(ells));
for n = 1:numel(ells)
    columns = (ell_of_state == ells(n));
    P{n} = V(:, columns) * V(:, columns)';
    D(n) = sum(columns);
end
end
