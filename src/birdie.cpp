#include "birdie.h"

// MAP estimate of Dirichlet-Multinomial posterior w/ columnwise prior alpha
// [[Rcpp::export]]
MatrixXd dirichlet_map(const IntegerVector Y, const MatrixXd r_probs,
                        const VectorXd prior_alpha) {
    int n_y = prior_alpha.size();
    int n_r = r_probs.cols();
    int N = r_probs.rows();

    MatrixXd post(n_y, n_r);
    for (int j = 0; j < n_r; j++) post.col(j) = prior_alpha.array() - 1.0; // -1 to help with the mode later

    for (int i = 0; i < N; i++) {
        post.row(Y[i] - 1) += r_probs.row(i);
    }

    post = post.array().rowwise() / post.colwise().sum().array();
    return post;
}

// [[Rcpp::export]]
List em_nocov(const IntegerVector Y, const Eigen::MatrixXd p_rxs,
              const Eigen::VectorXd prior_alpha, int iter) {
    MatrixXd est = dirichlet_map(Y, p_rxs, prior_alpha);
    MatrixXd est0 = est;

    MatrixXd p_ryxs(p_rxs.rows(), p_rxs.cols());

    for (int i = 0; i < iter; i++) {
        // E step
        p_ryxs = calc_bayes(Y, est, p_rxs);

        // M step
        est = dirichlet_map(Y, p_ryxs, prior_alpha);
    }

    return List::create(
        _["map"] = est,
        _["map0"] = est0,
        _["p_ryxs"] = p_ryxs
    );
}
