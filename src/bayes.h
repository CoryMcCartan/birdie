#ifndef BAYES_H
#define BAYES_H

#include <vector>

#include "birdie_types.h"

// [[Rcpp::export]]
Eigen::MatrixXd calc_bayes(
        const Eigen::VectorXi Y,
        const Eigen::VectorXi X,
        const std::vector<Eigen::MatrixXd> lik,
        const Eigen::MatrixXd prior,
        int n_x
);

std::vector<Eigen::MatrixXd> wrap_lik(const Eigen::MatrixXd lik);

#endif
