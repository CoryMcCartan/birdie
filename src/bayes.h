#ifndef BAYES_H
#define BAYES_H

#include <vector>

#include "birdie_types.h"

// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd calc_bayes(
        const Eigen::VectorXi Y,
        const Eigen::VectorXi X,
        const Eigen::VectorXd lik,
        const Eigen::MatrixXd prior,
        int n_x, int n_y
);

#endif
