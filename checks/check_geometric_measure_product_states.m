% CHECK_GEOMETRIC_MEASURE_PRODUCT_STATES  Inner maximum of Eq. (36) over product states: Neel-type
% states are optimal (Sec. 6, App. H).
% No CVX needed. 5-10 min. Expected: excess of the numerical search over the Neel value <= 1e-9.

%% Setup
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));

settings.random_starts = 10;
settings.restarts = 2;
settings.max_evaluations = 3000;
settings.seed = 0;
lambdas = -linspace(0.05, 1.5, 30);
test_values_w = [-0.2, -0.5, -0.85];

for N = [4 6]
    tic;
    [ells, ~, ~, ell_of_state, V] = total_spin_projectors(N);
    w = 2*ells.*(ells + 1)/N - 1;
    [~, p_neel] = neel_total_spin_distribution(N);

    %% 1. three ways to compute the largest eigenvalue
    rng(3);
    x = 3 * rand(1, 2*N);
    phi = product_state(x(1:N), x(N+1:2*N));
    lambda = -0.7;
    p = total_spin_weights(phi, V, ell_of_state, ells);
    small_matrix = largest_secular_root(p, lambda * w);
    secular = fzero(@(eta) sum(p ./ (eta - lambda*w)) - 1, [max(lambda*w) + 1e-9, max(lambda*w) + 10]);
    W_full = V * diag(2*ell_of_state.*(ell_of_state + 1)/N - 1) * V';
    direct = max(eig((lambda*W_full + phi*phi') + (lambda*W_full + phi*phi')') / 2);
    fprintf('\nN = %d: largest eigenvalue from the weights %.12f, secular equation %.12f, full matrix %.12f\n', ...
            N, small_matrix, secular, direct);

    %% 2. numerical maximization over product states, compared with the Neel state
    g_search = zeros(size(lambdas));
    g_neel = zeros(size(lambdas));
    for n = 1:numel(lambdas)
        g_search(n) = max_over_product_states(lambdas(n), N, V, ell_of_state, ells, settings);
        g_neel(n) = largest_secular_root(p_neel, lambdas(n) * w);
    end
    fprintf('N = %d: largest excess of the numerical search over the Neel value: %.2e\n', ...
            N, max(g_search - g_neel));

    %% 3. Legendre bound  E_G >= max over lambda of [1 + lambda*w - g(lambda)]
    for w_test = test_values_w
        bounds = 1 + lambdas * w_test - g_search;
        [best, k] = max(bounds);
        fprintf('N = %d, w = %+.2f: geometric-measure bound %.4f at lambda = %.3f  (BSA bound: %.2f)\n', ...
                N, w_test, best, lambdas(k), -w_test);
    end
    fprintf('(time for N = %d: %.0f s)\n', N, toc);
end
