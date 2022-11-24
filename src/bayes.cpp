#include "bayes.h"

// [[Rcpp::export(rng=false)]]
NumericMatrix calc_bayes_bisg(const IntegerVector S, const IntegerVector GX,
                              const NumericMatrix p_sr, const NumericMatrix p_gxr,
                              const NumericVector p_r) {
    int n_r = p_r.size();
    int n = S.size();
    NumericMatrix out(n, n_r);
    NumericVector sums(n);
    // do Bayes
    for (int j = 0; j < n_r; j++) {
        for (int i = 0; i < n; i++) {
            out(i, j) = p_sr(S[i] - 1, j) *  p_gxr(GX[i] - 1, j) * p_r[j];
            if (j == 0) {
                sums[i] = out(i, j);
            } else {
                sums[i] += out(i, j);
            }
        }
    }
    // normalize
    for (int i = 0; i < n_r; i++) {
        out(_, i) = out(_, i) / sums;
    }

    return out;
}

// exported in header
Eigen::MatrixXd calc_bayes(const Eigen::VectorXi Y, const Eigen::VectorXi X,
                           const Eigen::VectorXd lik,
                           const Eigen::MatrixXd prior, int n_x, int n_y) {
    int n_r = prior.cols();
    int N = Y.size();
    MatrixXd out(N, n_r);
    ArrayXd sums(N);
    // do Bayes
    for (int j = 0; j < n_r; j++) {
        for (int i = 0; i < N; i++) {
            out(i, j) = lik[est_idx(j, Y[i] - 1, X[i] - 1, n_r, n_y)] * prior(i, j);
            if (j == 0) {
                sums[i] = out(i, j);
            } else {
                sums[i] += out(i, j);
            }
        }
    }
    // normalize
    for (int i = 0; i < n_r; i++) {
        out.col(i) = out.col(i).array() / sums;
    }

    return out;
}
