function rho = gibbs_state(H, T)
% GIBBS_STATE  Thermal state exp(-H/T)/Z, Eq. (2).

[V, E] = eig(full(H));
E = diag(E);

if T == 0
    w = zeros(size(E));
    w(1) = 1;
else
    w = exp(-(E - E(1)) / T);
    w = w / sum(w);
end

rho = V * diag(w) * V';
rho = (rho + rho') / 2;
end
