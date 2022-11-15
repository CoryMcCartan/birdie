#include "bayes.h"

// [[Rcpp::export]]
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

// [[Rcpp::export]]
NumericMatrix calc_bayes(const IntegerVector Y,
                         const NumericMatrix lik, const NumericMatrix prior) {
    int n_r = prior.cols();
    int n = Y.size();
    NumericMatrix out(n, n_r);
    NumericVector sums(n);
    // do Bayes
    for (int j = 0; j < n_r; j++) {
        for (int i = 0; i < n; i++) {
            out(i, j) = lik(Y[i] - 1, j) *  prior(i, j);
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
