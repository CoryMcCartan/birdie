#ifndef TYPES_H
#define TYPES_H

// debug
#define PRINTLN Rcout << __LINE__ << " [" << __func__ << "(), " << __FILE__ << "]\n";

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;
using namespace arma;

#endif
