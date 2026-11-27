function value = max_over_ppt_weights(r, C)
% MAX_OVER_PPT_WEIGHTS  Linear program over the PPT total-spin weights (App. H).

require_cvx();
scale = max(abs(r));
r = r(:) / scale;
L = numel(r);

cvx_begin quiet
    variable t(L) nonnegative
    maximize( r' * t )
    subject to
        sum(t) == 1;
        C * t >= 0;
cvx_end

value = cvx_optval * scale;
end
