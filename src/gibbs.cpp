#include "gibbs.h"

// [[Rcpp::export]]
umat gibbs_race(int iter, int thin,
                const uvec X, const uvec S, const uvec GW,
                const mat lp_sr, const mat lp_wgr, const vec lp_r,
                const mat alpha) {
    // setup sizes and inits
    int N = X.size();
    int n_out = iter / thin;
    int n_r = lp_sr.n_cols;

    umat out(N, n_out);

    mat pr_baseline = calc_baseline_prob(n_r, N, S, GW, lp_sr, lp_wgr, lp_r);

    // initialize R
    uvec R(N);

    vec u = as<vec>(runif(N));
    for (int i = 0; i < N; i++)
        R[i] = rcatp(pr_baseline.col(i), u[i]);

    // run the Gibbs sampler
    for (int idx = 0; idx < iter; idx++) {
        if (idx % thin == 0) {
            out.col(idx/thin) = R;
        }
    }

    return out;
}


mat calc_baseline_prob(int n_r, int N, const uvec S, const uvec GW,
                       const mat lp_sr, const mat lp_wgr, const vec lp_r) {
    mat pr_baseline(n_r, N);
    for (int j = 0; j < n_r; j++) {
        for (int i = 0; i < N; i++) {
            pr_baseline(j, i) = lp_sr(S[i] - 1, j) + lp_wgr(GW[i] - 1, j) + lp_r[j];
        }
    }
    pr_baseline.each_row() -= max(pr_baseline, 0);
    pr_baseline = exp(pr_baseline);
    return pr_baseline;
}
