function e = energy_per_site_tdl(gamma, h, T)
% ENERGY_PER_SITE_TDL  Thermal energy per site of the infinite chain, Eq. (49).

G = majorana_G_tdl(3, gamma, h, T, 20000);
cxx = G(1 + 4);
cyy = G(-1 + 4);
mz  = -G(0 + 4);
e = -((1 + gamma)/2 * cxx + (1 - gamma)/2 * cyy) - h * mz;
end
