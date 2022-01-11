/***********************************
 * GIBBS SAMPLER FOR RACE IMPUTATION
 * Cory McCartan, Janurary 2022
 ***********************************/

#include "gibbs.h"

// [[Rcpp::export]]
umat gibbs_race(int iter, int thin,
                const uvec &X, const uvec &S, const uvec &GZ,
                const mat &lp_sr, const mat &lp_wgr, const vec &lp_r,
                const mat &alpha,
                int verbosity) {
    // setup sizes and inits
    int N = X.size();
    int n_out = iter / thin;
    iter -= (iter - 1) % thin;
    int n_r = lp_sr.n_cols;
    int n_gz = lp_wgr.n_rows;

    umat out(N, n_out);
    // initialize R
    mat pr_baseline = calc_baseline_prob(n_r, N, S, GZ, lp_sr, lp_wgr, lp_r);
    uvec R(N);
    vec alpha_p_r = 100 * exp(lp_r);
    vec p_r = exp(lp_r);

    vec u = as<vec>(runif(N));
    for (int i = 0; i < N; i++)
        R[i] = rcatp(pr_baseline.col(i) % p_r, u[i]);
    out.col(0) = R;

    // initialize counts
    vec alpha_r = sum(alpha, 1); // row sums: alpha is n_r by n_x
    lookup_GZRX n_gzrx;
    mat m_gzr(n_r, n_gz);
    for (int i = 0; i < n_gz; i++) m_gzr.col(i) = alpha_r;
    for (int i = 0; i < N; i++) {
        setup_n(n_gzrx, alpha, GZ[i], R[i], X[i], n_r);
        m_gzr(R[i] - 1, GZ[i] - 1)++;
    }
    std::vector<lookup_GZRX_it> it_gzx(N);
    for (int i = 0; i < N; i++) {
        const idx_GZRX i_n {GZ[i], X[i]};
        it_gzx[i] = n_gzrx.find(i_n);
    }
    // log cell counts to gauge sparsity
    if (verbosity >= 1) {
        Rcout << "GZ cells: " << m_gzr.size() << "\n";
        Rcout << "GZ/X cells: " << n_gzrx.size() << "\n";
    }

    // run the Gibbs sampler
    RObject bar = cli_progress_bar(iter, NULL);
    try {
        for (int it = 1; it < iter; it++) {
            Rcpp::checkUserInterrupt();
            if (CLI_SHOULD_TICK) cli_progress_set(bar, it);
            if (it % thin == 0) {
                out.col(it / thin) = R;
            }

            // R_i | R_{-i}
            vec u = as<vec>(runif(N));
            for (int i = 0; i < N; i++) {
                // compute n^{-i} and m^{-i}
                vec n_i = it_gzx[i]->second;
                vec m_i = m_gzr.col(GZ[i] - 1);
                n_i[R[i] - 1]--;
                m_i[R[i] - 1]--;
                // compute Pr(R_i | R_{-i}, G, S, X, W)
                vec probs = n_i / m_i % (pr_baseline.col(i) % p_r);
                // sample R_i
                R[i] = rcatp(probs, u[i]);
                // update n and m
                n_i[R[i] - 1]++;
                m_i[R[i] - 1]++;
            }

            // p_r | R
            p_r = rdirichlet(1, sum(m_gzr, 1) + alpha_p_r).row(0).t();
        }
    } catch (Rcpp::internal::InterruptedException e) {
        Rcerr << "Interrupted. Only partial results available.\n";
        cli_progress_done(bar);
        return out;
    }
    cli_progress_done(bar);

    return out;
}

void setup_n(lookup_GZRX &n_gzrx, const mat &alpha,
             int gz, int r, int x, int n_r) {
    const idx_GZRX i_n {gz, x};
    auto it_n = n_gzrx.find(i_n);

    if (it_n == n_gzrx.end()) {
        vec init = alpha.col(x - 1);
        init[r - 1]++;
        n_gzrx.insert(it_n, {i_n, init});
    } else {
        it_n->second[r - 1]++;
    }
}

mat calc_baseline_prob(int n_r, int N, const uvec &S, const uvec &GZ,
                       const mat &lp_sr, const mat &lp_wgr, const vec &lp_r) {
    mat pr_baseline(n_r, N);
    for (int j = 0; j < n_r; j++) {
        for (int i = 0; i < N; i++) {
            //pr_baseline(j, i) = lp_sr(S[i] - 1, j) + lp_wgr(GZ[i] - 1, j) + lp_r[j];
            pr_baseline(j, i) = lp_sr(S[i] - 1, j) + lp_wgr(GZ[i] - 1, j);
        }
    }
    pr_baseline.each_row() -= max(pr_baseline, 0); // column max
    pr_baseline = exp(pr_baseline);
    return pr_baseline;
}
