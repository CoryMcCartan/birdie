#ifndef TYPES_H
#define TYPES_H

// debug
#define PRINTLN Rcout << __LINE__ << " [" << __func__ << "(), " << __FILE__ << "]\n";

#include <RcppArmadillo.h>
#include <tuple>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;
using namespace arma;

typedef std::tuple<int, int> idx_GZRX;
typedef std::map<idx_GZRX, vec> lookup_GZRX;
typedef std::map<idx_GZRX, vec>::iterator lookup_GZRX_it;

#endif
