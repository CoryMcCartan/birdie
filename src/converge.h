#ifndef CONVERGE_H
#define CONVERGE_H

#include "birdie_types.h"

// [[Rcpp::export]]
bool check_convergence(
        Eigen::ArrayXd est,
        Eigen::ArrayXd last_est,
        double abstol,
        double reltol
);


#endif
