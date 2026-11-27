function Mhat = majorana_state_symbol_tdl(kind, varargin)
% MAJORANA_STATE_SYMBOL_TDL  Majorana correlations of the infinite chain in momentum space:
% thermal state, state at time t, steady state (GGE).

switch lower(kind)
    case 'thermal'
        [gamma, h, T] = varargin{1:3};
        k = optional_k(varargin, 4);
        X = 1i * majorana_symbol_tdl(gamma, h, k);
        if T > 0
            F = apply_function_2x2_hermitian(X, @(w) tanh(w / (2*T)));
        else
            F = apply_function_2x2_hermitian(X, @(w) sign(w));
        end
        Mhat = 1i * F;
    case 'evolved'
        [Mhat0, gamma_f, h_f, t] = varargin{1:4};
        k = optional_k(varargin, 5);
        X = 1i * majorana_symbol_tdl(gamma_f, h_f, k);
        U = apply_function_2x2_hermitian(X, @(w) exp(-1i * w * t));
        Mhat = zeros(size(Mhat0));
        for n = 1:size(Mhat0, 3)
            Mhat(:, :, n) = U(:, :, n) * Mhat0(:, :, n) * U(:, :, n)';
        end
    case 'steady'
        [Mhat0, gamma_f, h_f] = varargin{1:3};
        k = optional_k(varargin, 4);
        X = 1i * majorana_symbol_tdl(gamma_f, h_f, k);
        [~, Pp, Pm, degenerate] = apply_function_2x2_hermitian(X, @(w) w);
        Mhat = zeros(size(Mhat0));
        for n = 1:size(Mhat0, 3)
            if degenerate(n)
                Mhat(:, :, n) = Mhat0(:, :, n);
            else
                Mhat(:, :, n) = Pp(:, :, n) * Mhat0(:, :, n) * Pp(:, :, n) ...
                              + Pm(:, :, n) * Mhat0(:, :, n) * Pm(:, :, n);
            end
        end
    otherwise
        error('majorana_state_symbol_tdl: unknown kind "%s"', kind);
end
end

function k = optional_k(args, position)
if numel(args) >= position && ~isempty(args{position})
    k = args{position};
else
    k = momentum_grid_tdl();
end
end
