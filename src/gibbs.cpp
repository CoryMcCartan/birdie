/***********************************
 * GIBBS SAMPLER FOR RACE IMPUTATION
 * Cory McCartan, Janurary 2022
 ***********************************/

#include "gibbs.h"

// [[Rcpp::export]]
mat gibbs_race(int iter, int warmup,
               const uvec &X, const uvec &S, const uvec &GZ,
               const mat &M_sr, const mat &N_gzr,
               const mat &alpha_sr, const mat &beta_gzr,
               int verbosity) {
    // setup sizes and inits
    int N = X.size();
    int n_r = M_sr.n_cols;
    int n_s = M_sr.n_rows;
    int n_gz = N_gzr.n_rows;

    mat out(N, n_r, fill::zeros);

    // initialize R
    uvec R(N);
    vec p_r(n_r, fill::value(1.0 / n_r));
    vec u = as<vec>(runif(N));
    for (int i = 0; i < N; i++)
        R[i] = rcatp(p_r, u[i]);

    // initialize counts
    mat m_sr = (M_sr + alpha_sr).t();
    vec m_r = sum(m_sr, 1);
    mat n_gzr = (N_gzr + beta_gzr).t();
    for (int i = 0; i < N; i++) {
        m_sr(R[i] - 1, S[i] - 1)++;
        n_gzr(R[i] - 1, GZ[i] - 1)++;
        m_r[R[i] - 1]++;
    }

    // run the Gibbs sampler
    RObject bar = cli_progress_bar(iter, NULL);
    try {
        for (int it = 0; it < iter; it++) {
            Rcpp::checkUserInterrupt();
            if (CLI_SHOULD_TICK) cli_progress_set(bar, it);

            // R_i | R_{-i}
            vec u = as<vec>(runif(N));
            for (int i = 0; i < N; i++) {
                // compute n^{-i} and m^{-i}
                vec m_i = m_sr.col(S[i] - 1);
                vec n_i = n_gzr.col(GZ[i] - 1);
                n_i[R[i] - 1]--;
                m_i[R[i] - 1]--;
                m_r[R[i] - 1]--;
                // compute Pr(R_i | R_{-i}, G, S, Z)
                vec probs = n_i % m_i / m_r;

                // sample R_i
                R[i] = rcatp(probs, u[i]);
                // update n and m
                n_i[R[i] - 1]++;
                m_i[R[i] - 1]++;
                m_r[R[i] - 1]++;

                if (it >= warmup) {
                    out(i, R[i] - 1) += 1.0 / (iter - warmup);
                }
            }

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
                       const mat &lp_sr, const mat &lp_gzr, const vec &lp_r) {
    mat pr_baseline(n_r, N);
    for (int j = 0; j < n_r; j++) {
        for (int i = 0; i < N; i++) {
            pr_baseline(j, i) = lp_sr(S[i] - 1, j) + lp_gzr(GZ[i] - 1, j) + lp_r[j];
        }
    }
    pr_baseline.each_row() -= max(pr_baseline, 0); // column max
    pr_baseline = exp(pr_baseline);
    return pr_baseline;
}
