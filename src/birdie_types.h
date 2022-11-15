#ifndef TYPES_H
#define TYPES_H

// debug
#define PRINTLN Rcout << __LINE__ << " [" << __func__ << "(), " << __FILE__ << "]\n";
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;
using namespace Eigen;

#endif
