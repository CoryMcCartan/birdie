#ifndef TYPES_H
#define TYPES_H

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#define PRINTLN Rcout << __LINE__ << " [" << __func__ << "(), " << __FILE__ << "]\n";

#endif
