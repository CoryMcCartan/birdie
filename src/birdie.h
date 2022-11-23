#ifndef BIRDIE_H
#define BIRDIE_H

#include <vector>

#include "birdie_types.h"
#include "converge.h"
#include "bayes.h"

// [[Rcpp::export]]
std::vector<MatrixXd> dirichlet_map(
        const Eigen::VectorXi Y,
        const Eigen::VectorXi X,
        const Eigen::MatrixXd r_probs,
        const Eigen::VectorXd prior_alpha,
        int n_x
);


#endif
