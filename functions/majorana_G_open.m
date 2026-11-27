function G = majorana_G_open(N, gamma, h, T, J)
% MAJORANA_G_OPEN  Fermionic correlation matrix G_mn of an open chain in equilibrium (App. A).

if nargin < 5
    J = 1;
end

hm = zeros(N);
D  = zeros(N);
for j = 1:N-1
    hm(j, j+1) = J;
    hm(j+1, j) = J;
    D(j, j+1) =  J * gamma;
    D(j+1, j) = -J * gamma;
end
for j = 1:N
    hm(j, j) = -2 * h;
end

M = [hm, D; -D, -hm];
[W, e] = eig(M);
e = diag(e);
f = 1 ./ (1 + exp(-e / T));
Cm = W * diag(f) * W';

C11 = Cm(1:N, 1:N);        C12 = Cm(1:N, N+1:2*N);
C21 = Cm(N+1:2*N, 1:N);    C22 = Cm(N+1:2*N, N+1:2*N);
G = -(C21 + C22 - C11 - C12);
end
