// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "birdie_types.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// calc_bayes_bisg
NumericMatrix calc_bayes_bisg(const IntegerVector S, const IntegerVector GX, const NumericMatrix p_sr, const NumericMatrix p_gxr, const NumericVector p_r);
RcppExport SEXP _birdie_calc_bayes_bisg(SEXP SSEXP, SEXP GXSEXP, SEXP p_srSEXP, SEXP p_gxrSEXP, SEXP p_rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type S(SSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type GX(GXSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type p_sr(p_srSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type p_gxr(p_gxrSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type p_r(p_rSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_bayes_bisg(S, GX, p_sr, p_gxr, p_r));
    return rcpp_result_gen;
END_RCPP
}
// calc_bayes
Eigen::MatrixXd calc_bayes(const IntegerVector Y, const Eigen::MatrixXd lik, const Eigen::MatrixXd prior);
RcppExport SEXP _birdie_calc_bayes(SEXP YSEXP, SEXP likSEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type lik(likSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type prior(priorSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_bayes(Y, lik, prior));
    return rcpp_result_gen;
END_RCPP
}
// dirichlet_map
MatrixXd dirichlet_map(const IntegerVector Y, const MatrixXd r_probs, const VectorXd prior_alpha);
RcppExport SEXP _birdie_dirichlet_map(SEXP YSEXP, SEXP r_probsSEXP, SEXP prior_alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const MatrixXd >::type r_probs(r_probsSEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type prior_alpha(prior_alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(dirichlet_map(Y, r_probs, prior_alpha));
    return rcpp_result_gen;
END_RCPP
}
// em_nocov
List em_nocov(const IntegerVector Y, const Eigen::MatrixXd p_rxs, const Eigen::VectorXd prior_alpha, int iter);
RcppExport SEXP _birdie_em_nocov(SEXP YSEXP, SEXP p_rxsSEXP, SEXP prior_alphaSEXP, SEXP iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type p_rxs(p_rxsSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type prior_alpha(prior_alphaSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    rcpp_result_gen = Rcpp::wrap(em_nocov(Y, p_rxs, prior_alpha, iter));
    return rcpp_result_gen;
END_RCPP
}
// sum_grp
NumericVector sum_grp(const NumericVector x, const IntegerVector grp, int ngrp);
RcppExport SEXP _birdie_sum_grp(SEXP xSEXP, SEXP grpSEXP, SEXP ngrpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type grp(grpSEXP);
    Rcpp::traits::input_parameter< int >::type ngrp(ngrpSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_grp(x, grp, ngrp));
    return rcpp_result_gen;
END_RCPP
}
// sum_multi_grp
NumericMatrix sum_multi_grp(const IntegerVector x, const IntegerVector grp, const NumericVector wt, const NumericVector init, int nx, int ngrp);
RcppExport SEXP _birdie_sum_multi_grp(SEXP xSEXP, SEXP grpSEXP, SEXP wtSEXP, SEXP initSEXP, SEXP nxSEXP, SEXP ngrpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type grp(grpSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type init(initSEXP);
    Rcpp::traits::input_parameter< int >::type nx(nxSEXP);
    Rcpp::traits::input_parameter< int >::type ngrp(ngrpSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_multi_grp(x, grp, wt, init, nx, ngrp));
    return rcpp_result_gen;
END_RCPP
}
// gibbs_me
Eigen::MatrixXd gibbs_me(int iter, int warmup, const Eigen::VectorXi& S, const Eigen::VectorXi& GZ, const Eigen::MatrixXd& M_sr, const Eigen::MatrixXd& N_gzr, const Eigen::MatrixXd& alpha_gzr, const Eigen::MatrixXd& beta_sr, int verbosity);
RcppExport SEXP _birdie_gibbs_me(SEXP iterSEXP, SEXP warmupSEXP, SEXP SSEXP, SEXP GZSEXP, SEXP M_srSEXP, SEXP N_gzrSEXP, SEXP alpha_gzrSEXP, SEXP beta_srSEXP, SEXP verbositySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< int >::type warmup(warmupSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXi& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXi& >::type GZ(GZSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type M_sr(M_srSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type N_gzr(N_gzrSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type alpha_gzr(alpha_gzrSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type beta_sr(beta_srSEXP);
    Rcpp::traits::input_parameter< int >::type verbosity(verbositySEXP);
    rcpp_result_gen = Rcpp::wrap(gibbs_me(iter, warmup, S, GZ, M_sr, N_gzr, alpha_gzr, beta_sr, verbosity));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_birdie_calc_bayes_bisg", (DL_FUNC) &_birdie_calc_bayes_bisg, 5},
    {"_birdie_calc_bayes", (DL_FUNC) &_birdie_calc_bayes, 3},
    {"_birdie_dirichlet_map", (DL_FUNC) &_birdie_dirichlet_map, 3},
    {"_birdie_em_nocov", (DL_FUNC) &_birdie_em_nocov, 4},
    {"_birdie_sum_grp", (DL_FUNC) &_birdie_sum_grp, 3},
    {"_birdie_sum_multi_grp", (DL_FUNC) &_birdie_sum_multi_grp, 6},
    {"_birdie_gibbs_me", (DL_FUNC) &_birdie_gibbs_me, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_birdie(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
