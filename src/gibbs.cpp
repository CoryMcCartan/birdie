#include "gibbs.h"

// [[Rcpp::export]]
mat gibbs_me(int iter, int warmup, const uvec &S, const uvec &GZ,
             const mat &M_sr, const mat &N_gzr,
             const mat &alpha_gzr, const mat &beta_sr,
             int verbosity) {
    // setup sizes and inits
    int N = S.size();
    int n_r = M_sr.n_cols;
    // int n_s = M_sr.n_rows;
    // int n_gz = N_gzr.n_rows;

    mat out(N, n_r, fill::zeros);

    // initialize R
    uvec R(N);
    vec p_r(n_r, fill::value(1.0 / n_r));
    vec u = as<vec>(runif(N));
    for (int i = 0; i < N; i++)
        R[i] = rcatp(p_r, u[i]);

    // initialize counts
    mat m_sr = (M_sr + beta_sr).t();
    vec m_r = sum(m_sr, 1);
    mat n_gzr = (N_gzr + alpha_gzr).t();
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
