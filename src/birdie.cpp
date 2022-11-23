#include "birdie.h"

// EM model with discrete fixed effect
// Can fit no-pooling or complete pooling model depending on imputs
// [[Rcpp::export]]
List em_fixef(const Eigen::VectorXi Y, const Eigen::VectorXi X,
              const Eigen::MatrixXd p_rxs,
              const Eigen::VectorXd prior_alpha, int n_x, int iter,
              double abstol, double reltol) {
    int N = Y.size();
    int n_y = prior_alpha.size();
    std::vector<MatrixXd> ests = dirichlet_map(Y, X, p_rxs, prior_alpha, n_x);

    MatrixXd p_ryxs(p_rxs.rows(), p_rxs.cols());
    VectorXi ones = VectorXi::Constant(N, 1);
    VectorXd unif_prior = VectorXd::Constant(n_y, 1.0);

    MatrixXd est;
    MatrixXd last_est;

    int i;
    for (i = 0; i < iter; i++) {
        last_est = est;

        // E step
        p_ryxs = calc_bayes(Y, X, ests, p_rxs, n_x);

        // M step
        ests = dirichlet_map(Y, X, p_ryxs, prior_alpha, n_x);
        est = dirichlet_map(Y, ones, p_ryxs, unif_prior, 1)[0];

        if (i > 0 && check_convergence(est, last_est, abstol, reltol)) break;
    }

    // final E step
    p_ryxs = calc_bayes(Y, X, ests, p_rxs, n_x);

    return List::create(
        _["map"] = est,
        _["map_x"] = ests,
        _["p_ryxs"] = p_ryxs,
        _["iters"] = i + 1,
        _["converge"] = i < iter - 1
    );
}


// MAP estimate of Dirichlet-Multinomial posterior by X
// exported in header
std::vector<MatrixXd> dirichlet_map(
        const Eigen::VectorXi Y, const Eigen::VectorXi X,
        const Eigen::MatrixXd r_probs, const Eigen::VectorXd prior_alpha, int n_x) {
    int n_y = prior_alpha.size();
    int n_r = r_probs.cols();
    int N = r_probs.rows();

    std::vector<MatrixXd> post(n_x);
    for (int i = 0; i < n_x; i++) {
        post[i] = MatrixXd::Zero(n_y, n_r);
        for (int j = 0; j < n_r; j++) {
            post[i].col(j) = prior_alpha.array() - 1.0; // -1 to help with the mode later
        }
    }

    for (int i = 0; i < N; i++) {
        post[X[i] - 1].row(Y[i] - 1) += r_probs.row(i);
    }

    for (int i = 0; i < n_x; i++) {
        post[i] = post[i].array().rowwise() / post[i].colwise().sum().array();
    }

    return post;
}

// [[Rcpp::export]]
NumericMatrix sum_grp(const IntegerVector x, const IntegerVector grp,
                      const NumericVector wt, const NumericVector init,
                      int nx, int ngrp) {
    NumericMatrix out(ngrp, nx);

    for (int j = 0; j < ngrp; j++) {
        for (int k = 0; k < nx; k++) {
            out(j, k) = init[k];
        }
    }

    int N = x.size();
    for (int i = 0; i < N; i++) {
        out(grp[i] - 1, x[i] - 1) += wt[i];
    }

    return out;
}
