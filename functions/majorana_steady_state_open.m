function Mss = majorana_steady_state_open(M0, Af)
% MAJORANA_STEADY_STATE_OPEN  Steady state (GGE) of an open chain after a quench.

X = 1i * Af;
X = (X + X') / 2;
[V, w] = eig(X);
w = real(diag(w));
Mp = V' * M0 * V;
same = abs(w - w.') < 1e-9;
Mss = real(V * (Mp .* same) * V');
Mss = (Mss - Mss.') / 2;
end
