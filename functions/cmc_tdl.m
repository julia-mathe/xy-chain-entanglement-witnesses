function [omega, s, qmin, converged, R, curves] = cmc_tdl(gamma, h, T)
% CMC_TDL  omega_CMC of the infinite chain in equilibrium, Eqs. (25) and (29).

R = 80;
Rmax = 400;
tolerance = 1e-9;

while true
    [mz, Cxx, Cyy, Czz] = spin_correlations_tdl(R, gamma, h, T);
    Czc = Czz - mz^2;
    Czc(1) = 1 - mz^2;
    tail = max(abs([Cxx(end), Cyy(end), Czc(end)]));
    if tail < tolerance || R >= Rmax
        break
    end
    R = min(2*R, Rmax);
end
converged = tail < tolerance;

q = linspace(0, pi, 721);
C = {Cxx, Cyy, Czc};
S = zeros(3, numel(q));
Spi = zeros(1, 3);
S0 = zeros(1, 3);
for alpha = 1:3
    S(alpha, :) = structure_factor_from_correlations(C{alpha}, q);
    Spi(alpha) = staggered_structure_factor_euler(C{alpha});
    S0(alpha) = C{alpha}(1) + 2 * sum(C{alpha}(2:end));
end

s = zeros(1, 3);
qmin = zeros(1, 3);
if converged
    [s, k] = min(S, [], 2);
    s = s';
    qmin = q(k);
else
    for alpha = 1:3
        if alpha == 1
            candidates = [Inf, Spi(1)];
        else
            candidates = [S0(alpha), Spi(alpha)];
        end
        [s(alpha), k] = min(candidates);
        qmin(alpha) = (k - 1) * pi;
    end
end

omega = theorem1_omega(s);

curves.q = q;
curves.S = S;
curves.Spi = Spi;
curves.S0 = S0;
curves.mz = mz;
curves.tail = tail;
end
