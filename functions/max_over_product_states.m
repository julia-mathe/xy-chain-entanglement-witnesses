function [g, best_theta, best_phi] = max_over_product_states(lambda, N, V, ell_of_state, ells, settings)
% MAX_OVER_PRODUCT_STATES  Numerical search for the inner maximum of Eq. (36) over product states.

w = 2*ells.*(ells + 1)/N - 1;
a = lambda * w;

minus_objective = @(x) -largest_secular_root( ...
    total_spin_weights(product_state(x(1:N), x(N+1:2*N)), V, ell_of_state, ells), a);

if isfield(settings, 'seed')
    rng(settings.seed);
end
starts = zeros(settings.random_starts + 1, 2*N);
starts(1, :) = [repmat([0, pi], 1, N/2), zeros(1, N)];
for n = 1:settings.random_starts
    starts(n + 1, :) = [acos(2*rand(1, N) - 1), 2*pi*rand(1, N)];
end

options = optimset('MaxFunEvals', settings.max_evaluations, 'MaxIter', settings.max_evaluations, ...
                   'TolX', 1e-9, 'TolFun', 1e-12, 'Display', 'off');
g = -Inf;
for n = 1:size(starts, 1)
    x = starts(n, :);
    for repeat = 1:settings.restarts + 1
        [x, value] = fminsearch(minus_objective, x, options);
    end
    if -value > g
        g = -value;
        best_theta = x(1:N);
        best_phi = x(N+1:2*N);
    end
end
end
