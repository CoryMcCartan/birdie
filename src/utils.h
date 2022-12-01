#ifndef UTILS_H
#define UTILS_H

#include <cmath>

#include "birdie_types.h"

// [[Rcpp::export(rng = false)]]
bool check_convergence(
        Eigen::ArrayXd est,
        Eigen::ArrayXd last_est,
        double abstol,
        double reltol
);

// [[Rcpp::export(rng = false)]]
Eigen::VectorXd to_simplex(Eigen::VectorXd y);

// [[Rcpp::export(rng = false)]]
Eigen::VectorXd from_simplex(Eigen::VectorXd x);


#endif
