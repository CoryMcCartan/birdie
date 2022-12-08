#include <cmath>
#include "birdie.h"

// MAP estimate of Dirichlet-Multinomial posterior by X
// Works on a 3D array represented as a long vector (see `birdie_types.h` for indexing fns)
// [[Rcpp::export(rng=false)]]
Eigen::VectorXd dirichlet_map(
        const Eigen::VectorXi Y, const Eigen::VectorXi X,
        const Eigen::MatrixXd p_rxs, const Eigen::MatrixXd prior_yr, int n_x) {
    int N = Y.rows();
    int n_y = prior_yr.rows();
    int n_r = p_rxs.cols();

    Eigen::VectorXd post(n_x * n_y * n_r);
    Eigen::VectorXd prior_sum = prior_yr.colwise().sum().array() - n_y;
    Eigen::MatrixXd sums(n_r, n_x);
    for (int i = 0; i < n_x; i++) {
        sums.col(i) = prior_sum;
        for (int j = 0; j < n_y; j++) {
            for (int k = 0; k < n_r; k++) {
                post[est_idx(k, j, i, n_r, n_y)] = prior_yr(j, k) - 1.0; // -1 to help with the mode later
            }
        }
    }

    for (int i = 0; i < N; i++) {
        sums.col(X[i] - 1) += p_rxs.row(i);
        post.segment(est_col(Y[i] - 1, X[i] - 1, n_r, n_y), n_r) += p_rxs.row(i);
    }

    // normalize
    for (int i = 0; i < n_x; i++) {
        for (int j = 0; j < n_y; j++) {
            post.segment(est_col(j, i, n_r, n_y), n_r).array() /= sums.col(i).array();
        }
    }

    return post;
}

// combined calc_bayes + dirichlet_map (doesn't store intermediate p_ryxs matrix)
// input (Y|R) tables in `curr`
// if `sum_only` is true, returns aggregated data rather than Dirichlet mode
// [[Rcpp::export(rng=false)]]
Eigen::VectorXd em_dirichlet(
        const Eigen::VectorXd curr,
        const Eigen::VectorXi Y, const Eigen::VectorXi X,
        const Eigen::MatrixXd p_rxs, const Eigen::MatrixXd prior_yr,
        int n_x, bool sum_only=false) {
    int N = Y.size();
    int n_y = prior_yr.rows();
    int n_r = p_rxs.cols();

    // initialize output with prior
    VectorXd post(n_x * n_y * n_r);
    MatrixXd sums(n_r, n_x);
    if (sum_only) {
        post = VectorXd::Zero(n_x * n_y * n_r);
        sums = MatrixXd::Zero(n_r, n_x);
    } else {
        VectorXd prior_sum = prior_yr.colwise().sum().array() - n_y;
        for (int i = 0; i < n_x; i++) {
            sums.col(i) = prior_sum;
            for (int j = 0; j < n_y; j++) {
                for (int k = 0; k < n_r; k++) {
                    post[est_idx(k, j, i, n_r, n_y)] = prior_yr(j, k) - 1; // -1 if we are finding mode/MAP
                }
            }
        }
    }

    // do Bayes and sum
    VectorXd pr_i(n_r);
    for (int i = 0; i < N; i++) {
        int idx = est_col(Y[i] - 1, X[i] - 1, n_r, n_y);
        for (int k = 0; k < n_r; k++) {
            pr_i[k] = curr[idx + k] * p_rxs(i, k);
        }
        pr_i /= pr_i.sum();

        post.segment(idx, n_r) += pr_i;
        sums.col(X[i] - 1) += pr_i;
    }

    // convert to MLE
    if (!sum_only) {
        for (int i = 0; i < n_x; i++) {
            for (int j = 0; j < n_y; j++) {
                post.segment(est_col(j, i, n_r, n_y), n_r).array() /= sums.col(i).array();
            }
        }
    }

    return post;
}

// weighted version
// [[Rcpp::export(rng=false)]]
Eigen::VectorXd em_dirichlet_wt(
        const Eigen::VectorXd curr,
        const Eigen::VectorXi Y, const Eigen::VectorXi X,
        const Eigen::VectorXd wt,
        const Eigen::MatrixXd p_rxs, const Eigen::MatrixXd prior_yr, int n_x) {
    int N = Y.size();
    int n_y = prior_yr.rows();
    int n_r = p_rxs.cols();

    // initialize output with prior
    VectorXd post(n_x * n_y * n_r);
    MatrixXd sums(n_r, n_x);
    VectorXd prior_sum = prior_yr.colwise().sum().array() - n_y;
    for (int i = 0; i < n_x; i++) {
        sums.col(i) = prior_sum;
        for (int j = 0; j < n_y; j++) {
            for (int k = 0; k < n_r; k++) {
                post[est_idx(k, j, i, n_r, n_y)] = prior_yr(j, k) - 1; // -1 if we are finding mode/MAP
            }
        }
    }

    // do Bayes and sum
    VectorXd pr_i(n_r);
    for (int i = 0; i < N; i++) {
        int idx = est_col(Y[i] - 1, X[i] - 1, n_r, n_y);
        for (int k = 0; k < n_r; k++) {
            pr_i[k] = curr[idx + k] * p_rxs(i, k);
        }
        pr_i /= pr_i.sum();

        post.segment(idx, n_r) += pr_i * wt[i];
        sums.col(X[i] - 1) += pr_i * wt[i];
    }

    for (int i = 0; i < n_x; i++) {
        for (int j = 0; j < n_y; j++) {
            post.segment(est_col(j, i, n_r, n_y), n_r).array() /= sums.col(i).array();
        }
    }

    return post;
}

// multiplication helper
// [[Rcpp::export(rng=false)]]
Eigen::VectorXd resid_mult(const Eigen::VectorXd m_coef, const Eigen::VectorXi idxs,
                           const Eigen::MatrixXd r_probs, int k, int n_k) {
    int N = r_probs.rows();
    int n_r = r_probs.cols();
    Eigen::VectorXd out(N);

    k--; // to zero-index
    for (int i = 0; i < N; i++) {
        out[i] = 0;
        for (int j = 0; j < n_r; j++) {
            int idx = k + n_k * (j + n_r * (idxs[i] - 1));
            out[i] += m_coef[idx] * r_probs(i, j);
        }
    }

    return out;
}
