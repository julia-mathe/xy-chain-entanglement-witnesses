function k = momentum_grid_tdl(nk)
% MOMENTUM_GRID_TDL  Momenta k = 2 pi (m + 1/2) / nk used for the infinite chain.

if nargin < 1
    nk = 4096;
end
k = 2*pi*((0:nk-1) + 0.5) / nk;
end
