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


#endif
