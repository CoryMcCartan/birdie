#include "gibbs.h"


// [[Rcpp::export]]
Eigen::MatrixXd gibbs_me(
        int iter, int warmup, const Eigen::VectorXi &S, const Eigen::VectorXi &GZ,
        const Eigen::MatrixXd &M_sr, const Eigen::MatrixXd &N_gzr,
        const Eigen::MatrixXd &alpha_gzr, const Eigen::MatrixXd &beta_sr,
        int cores=0, int verbosity=3) {
    // setup sizes and inits
    int N = S.size();
    int n_r = M_sr.cols();

    MatrixXd out = MatrixXd::Zero(N, n_r);

    // initialize R
    ArrayXi R(N);
    ArrayXd p_r = ArrayXd::Constant(n_r, 1.0 / n_r);
    VectorXd u = as<VectorXd>(runif(N));
    for (int i = 0; i < N; i++)
        R[i] = rcatp(p_r, u[i]);

    // initialize counts
    MatrixXd m_sr = (M_sr + beta_sr).transpose();
    ArrayXd m_r = m_sr.rowwise().sum();
    MatrixXd n_gzr = (N_gzr + alpha_gzr).transpose();
    for (int i = 0; i < N; i++) {
        m_sr(R[i] - 1, S[i] - 1)++;
        n_gzr(R[i] - 1, GZ[i] - 1)++;
        m_r[R[i] - 1]++;
    }

    // run the Gibbs sampler
    RObject bar = cli_progress_bar(iter, NULL);
    RcppThread::ThreadPool pool(cores);
    try {
        for (int it = 0; it < iter; it++) {
            Rcpp::checkUserInterrupt();
            if (CLI_SHOULD_TICK) cli_progress_set(bar, it);

            // R_i | R_{-i}
            VectorXd u = as<VectorXd>(runif(N));
            pool.parallelFor(0, N, [&] (int i) {
                // compute n^{-i} and m^{-i}
                ArrayXd m_i = m_sr.col(S[i] - 1);
                ArrayXd n_i = n_gzr.col(GZ[i] - 1);
                n_i[R[i] - 1]--;
                m_i[R[i] - 1]--;
                m_r[R[i] - 1]--;
                // compute Pr(R_i | R_{-i}, G, S, Z)
                ArrayXd probs = n_i * m_i / m_r;

                // sample R_i
                R[i] = rcatp(probs, u[i]);
                // update n and m
                n_i[R[i] - 1]++;
                m_i[R[i] - 1]++;
                m_r[R[i] - 1]++;

                if (it >= warmup) {
                    out(i, R[i] - 1) += 1.0 / (iter - warmup);
                }
            });
            pool.wait();
        }
    } catch (Rcpp::internal::InterruptedException e) {
        Rcerr << "Interrupted. Only partial results available.\n";
        cli_progress_done(bar);
        return out;
    }
    cli_progress_done(bar);

    return out;
}
