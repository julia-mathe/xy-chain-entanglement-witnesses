function v = ppt_negative_part(rho, A, N)
% PPT_NEGATIVE_PART  v_PPT = Tr(W_A rho): sum of the negative eigenvalues of rho^{T_A} (Sec. 2.2.1).
X = partial_transpose(rho, A, N);
d = real(eig(full((X + X') / 2)));
v = sum(d(d < -1e-12));
end
