function Gamma = covariance_window_open(G, first_site, L)
% COVARIANCE_WINDOW_OPEN  Covariance matrix of L consecutive spins of an open chain (App. A).

sites = first_site : first_site + L - 1;
mz = -diag(G(sites, sites));
Gamma = zeros(3*L);

for a = 1:L
    for b = 1:L
        i = sites(a);   j = sites(b);
        ia = 3*(a-1);   ib = 3*(b-1);
        if i == j
            Gamma(ia+1, ib+1) = 1;
            Gamma(ia+2, ib+2) = 1;
            Gamma(ia+3, ib+3) = 1 - mz(a)^2;
            continue
        end
        lo = min(i, j);   r = abs(i - j);
        Mx = zeros(r);
        My = zeros(r);
        for p = 1:r
            for k = 1:r
                Mx(p, k) = G(lo + p - 1, lo + k);
                My(p, k) = G(lo + p, lo + k - 1);
            end
        end
        Gamma(ia+1, ib+1) = det(Mx);
        Gamma(ia+2, ib+2) = det(My);
        Gamma(ia+3, ib+3) = G(i,i)*G(j,j) - G(i,j)*G(j,i) - mz(a)*mz(b);
    end
end
end
