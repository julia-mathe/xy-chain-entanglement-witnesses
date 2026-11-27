function map = rdbu_colormap(n)
% RDBU_COLORMAP  Diverging blue-white-red colour map.

if nargin < 1
    n = 256;
end
reference = [  5  48  97
              33 102 172
              67 147 195
             146 197 222
             209 229 240
             247 247 247
             253 219 199
             244 165 130
             214  96  77
             178  24  43
             103   0  31] / 255;
positions = linspace(0, 1, size(reference, 1));
map = interp1(positions, reference, linspace(0, 1, n));
end
