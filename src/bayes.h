#ifndef BAYES_H
#define BAYES_H

#include "birdie_types.h"

// [[Rcpp::export]]
Eigen::MatrixXd calc_bayes(const IntegerVector Y,
                           const Eigen::MatrixXd lik, const Eigen::MatrixXd prior);

#endif
