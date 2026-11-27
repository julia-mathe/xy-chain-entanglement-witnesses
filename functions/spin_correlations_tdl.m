function [mz, Cxx, Cyy, Czz] = spin_correlations_tdl(R, gamma, h, T, nk)
% SPIN_CORRELATIONS_TDL  Thermal correlations of the infinite chain, Eqs. (41)-(44).

if nargin < 5
    nk = [];
end

G  = majorana_G_tdl(R + 1, gamma, h, T, nk);
Gr = @(r) G(r + (R + 1) + 1);

mz = -Gr(0);

Czz = zeros(1, R + 1);
Czz(1) = 1;
for r = 1:R
    Czz(r + 1) = Gr(0)^2 - Gr(r) * Gr(-r);
end

[a, b] = ndgrid(1:R, 1:R);
Mx = Gr(b - a + 1);
My = Gr(b - a - 1);

Cxx = [1, leading_principal_minors(Mx)];
Cyy = [1, leading_principal_minors(My)];
end
