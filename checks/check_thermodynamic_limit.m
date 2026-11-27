% CHECK_THERMODYNAMIC_LIMIT  omega_CMC of the infinite chain from Eqs. (25) and (29) (Table 1, App. A).
% No CVX needed. Expected: XX point 0.18940 (s = 0.1574 0.1574 0.0639), Ising point 0.55870
% (s = 0.2986 0.4528 0.3660); Toeplitz differences ~1e-16; accelerated S_x(pi) independent of R.

%% Setup
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));

%% Part 1: omega_CMC of the infinite chain
points = [0    0    0.1
          1    1    0.1
          0.5  0.5  0.5
          0.5  0.5  1.0
          1    0.3  0.5
          0    1.5  1.0
          0    2    0.8];
fprintf('Part 1: optimal CMC value of the infinite chain\n');
fprintf('%6s %6s %6s %10s %8s %8s %8s   %-12s %10s %5s\n', 'gamma', 'h', 'T', 'omega', 's_x', 's_y', ...
        's_z', 'q/pi', 'converged', 'R');
for n = 1:size(points, 1)
    [omega, s, qmin, converged, R] = cmc_tdl(points(n,1), points(n,2), points(n,3));
    if converged
        answer = 'yes';
    else
        answer = 'no';
    end
    fprintf('%6.2f %6.2f %6.2f %10.5f %8.4f %8.4f %8.4f   %-12s %10s %5d\n', points(n,:), omega, s, ...
            sprintf('%g ', qmin / pi), answer, R);
end

%% Part 2: all Toeplitz determinants at once versus separate determinants
fprintf('\nPart 2: Toeplitz determinants (R = 120)\n');
for n = [2 1 3]
    gamma = points(n,1);  h = points(n,2);  T = points(n,3);
    R = 120;
    [~, Cxx, Cyy] = spin_correlations_tdl(R, gamma, h, T);
    G = majorana_G_tdl(R + 1, gamma, h, T);
    Gr = @(r) G(r + (R + 1) + 1);
    [a, b] = ndgrid(1:R, 1:R);
    Mx = Gr(b - a + 1);
    My = Gr(b - a - 1);
    direct_x = zeros(1, R);
    direct_y = zeros(1, R);
    for r = 1:R
        direct_x(r) = det(Mx(1:r, 1:r));
        direct_y(r) = det(My(1:r, 1:r));
    end
    fprintf('  gamma = %.1f, h = %.1f, T = %.1f:  largest difference  xx %.1e,  yy %.1e\n', gamma, h, T, ...
            max(abs(Cxx(2:end) - direct_x)), max(abs(Cyy(2:end) - direct_y)));
end

%% Part 3: accelerated alternating sum for the staggered structure factor S_x(pi)
fprintf('\nPart 3: S_x(pi) with a plain and with an accelerated sum\n');
examples = [1 0.5 0.3
            1 0.5 0.5
            1 0.3 0.4
            1 1.0 0.1];
for n = 1:size(examples, 1)
    for R = [100 200 400 800]
        [~, Cxx] = spin_correlations_tdl(R, examples(n,1), examples(n,2), examples(n,3));
        plain = Cxx(1) + 2 * sum((-1).^(1:R) .* Cxx(2:end));
        accelerated = staggered_structure_factor_euler(Cxx);
        fprintf('  gamma = %g, h = %g, T = %g, R = %3d:  plain %.6f   accelerated %.6f   C_xx(R) = %.2e\n', ...
                examples(n,:), R, plain, accelerated, Cxx(end));
    end
end
