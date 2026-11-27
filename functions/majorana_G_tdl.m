function G = majorana_G_tdl(rmax, gamma, h, T, nk)
% MAJORANA_G_TDL  G_r of the infinite chain, Eq. (45).

if nargin < 5 || isempty(nk)
    if T >= 0.05
        nk = 4000;
    else
        nk = 40000;
    end
end

k = ((1:nk) - 0.5) * pi / nk;
Lambda = sqrt((cos(k) - h).^2 + gamma^2 * sin(k).^2);
if T > 0
    thermal_factor = tanh(Lambda / T);
else
    thermal_factor = ones(size(k));
end
A = (cos(k) - h) ./ Lambda .* thermal_factor;
B = gamma * sin(k) ./ Lambda .* thermal_factor;

r = -rmax:rmax;
G = zeros(size(r));
for n = 1:numel(r)
    G(n) = mean( cos(r(n) * k) .* A + sin(r(n) * k) .* B );
end
end
