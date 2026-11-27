function M = majorana_thermal_state(A, T)
% MAJORANA_THERMAL_STATE  M = i tanh(i A / 2T) of the Gibbs state.

X = 1i * A;
X = (X + X') / 2;
[V, w] = eig(X);
w = real(diag(w));
if T > 0
    f = tanh(w / (2*T));
else
    f = sign(w);
end
M = 1i * (V * diag(f) * V');
M = real(M);
M = (M - M.') / 2;
end
