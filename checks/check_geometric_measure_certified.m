% CHECK_GEOMETRIC_MEASURE_CERTIFIED  Lower bounds on the geometric measure, Eq. (36) (Table 3, App. H).
% Needs CVX. A few minutes. Expected (N = 10): omega = 0.154: certified 0.3993, Neel 0.4009;
% omega = 0.583: certified 0.0302, Neel 0.0620.

%% Setup
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));

%% Part 1: bounds assuming that Neel states are optimal
fprintf('Part 1: bounds assuming Neel optimality\n');
for N = [4 6 8 10]
    fprintf('  N = %2d:  omega = 0.15371 -> E_G >= %.4f,   omega = 0.5827 -> E_G >= %.4f\n', N, ...
            geometric_bound_neel(0.15371, N), geometric_bound_neel(0.5827, N));
end
[bound_xx, c_xx] = geometric_bound_large_N(0.189);
[bound_is, c_is] = geometric_bound_large_N(0.559);
fprintf('  large N: omega = 0.189 -> E_G >= %.4f (c = %.4f),  omega = 0.559 -> E_G >= %.4f (c = %.4f)\n', ...
        bound_xx, c_xx, bound_is, c_is);
fprintf(['  (for comparison, the approximation (1-omega)^2 (2-omega)/4, valid only for omega\n' ...
         '   close to 1, gives %.4f and %.4f)\n'], (1-0.189)^2*(2-0.189)/4, (1-0.559)^2*(2-0.559)/4);

%% Part 2: certified bounds (PPT relaxation)
fprintf('\nPart 2: certified bounds from the PPT relaxation\n');
cases = { 6,  0.154
          8,  0.154
          10, [0.154, 0.583] };
for n = 1:size(cases, 1)
    N = cases{n, 1};
    tic;
    [C, ells, D, report] = ppt_weight_constraints(N);
    fprintf('  N = %d: PPT conditions computed in %.0f s\n', N, toc);
    for k = 1:size(report, 1)
        fprintf('     cut size k = %d: largest off-diagonal element %.1e (should be ~0), %d distinct conditions\n', ...
                k, report(k, 1), report(k, 2));
    end
    for omega = cases{n, 2}
        [bound, lambda_best, eta_best] = geometric_bound_certified(omega, N, C, ells);
        fprintf('     omega = %.3f: certified E_G >= %.4f at lambda = %.3f (g_up = %.6f);  Neel value %.4f\n', ...
                omega, bound, lambda_best, eta_best, geometric_bound_neel(omega, N));

        %% Part 3: cross-check of g_up(lambda) by linear programming + bisection in eta
        w = 2*ells.*(ells + 1)/N - 1;
        a = lambda_best * w;
        excess = @(eta) max_over_ppt_weights(1 ./ (eta - a), C) - 1;
        low = max(a) + 1e-9;
        high = max(a) + 1;
        while excess(high) > 0
            high = max(a) + 2 * (high - max(a));
        end
        g_up = fzero(excess, [low, high], optimset('TolX', 1e-10));
        fprintf('                   cross-check: g_up from LP + bisection = %.6f\n', g_up);
    end
end
