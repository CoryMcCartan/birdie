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

// indexing functions for 3D estimate arrays
inline int est_idx(int r, int y, int grp, int n_r, int n_y) {
    return r + n_r*(y + n_y*grp);
}
inline int est_col(int y, int grp, int n_r, int n_y) {
    return n_r*(y + n_y*grp);
}


#endif
