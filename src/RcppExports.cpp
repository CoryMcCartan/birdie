// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "birdie_types.h"
#include <RcppEigen.h>
#include <RcppThread.h>
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
Eigen::MatrixXd calc_bayes(const Eigen::VectorXi Y, const Eigen::VectorXi X, const Eigen::VectorXd lik, const Eigen::MatrixXd prior, int n_x, int n_y);
RcppExport SEXP _birdie_calc_bayes(SEXP YSEXP, SEXP XSEXP, SEXP likSEXP, SEXP priorSEXP, SEXP n_xSEXP, SEXP n_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXi >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXi >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type lik(likSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< int >::type n_x(n_xSEXP);
    Rcpp::traits::input_parameter< int >::type n_y(n_ySEXP);
    rcpp_result_gen = Rcpp::wrap(calc_bayes(Y, X, lik, prior, n_x, n_y));
    return rcpp_result_gen;
END_RCPP
}
// dirichlet_map
Eigen::VectorXd dirichlet_map(const Eigen::VectorXi Y, const Eigen::VectorXi X, const Eigen::MatrixXd p_rxs, const Eigen::MatrixXd prior_yr, int n_x);
RcppExport SEXP _birdie_dirichlet_map(SEXP YSEXP, SEXP XSEXP, SEXP p_rxsSEXP, SEXP prior_yrSEXP, SEXP n_xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXi >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXi >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type p_rxs(p_rxsSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type prior_yr(prior_yrSEXP);
    Rcpp::traits::input_parameter< int >::type n_x(n_xSEXP);
    rcpp_result_gen = Rcpp::wrap(dirichlet_map(Y, X, p_rxs, prior_yr, n_x));
    return rcpp_result_gen;
END_RCPP
}
// em_dirichlet
Eigen::VectorXd em_dirichlet(const Eigen::VectorXd curr, const Eigen::VectorXi Y, const Eigen::VectorXi X, const Eigen::MatrixXd p_rxs, const Eigen::MatrixXd prior_yr, int n_x, bool sum_only);
RcppExport SEXP _birdie_em_dirichlet(SEXP currSEXP, SEXP YSEXP, SEXP XSEXP, SEXP p_rxsSEXP, SEXP prior_yrSEXP, SEXP n_xSEXP, SEXP sum_onlySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type curr(currSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXi >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXi >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type p_rxs(p_rxsSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type prior_yr(prior_yrSEXP);
    Rcpp::traits::input_parameter< int >::type n_x(n_xSEXP);
    Rcpp::traits::input_parameter< bool >::type sum_only(sum_onlySEXP);
    rcpp_result_gen = Rcpp::wrap(em_dirichlet(curr, Y, X, p_rxs, prior_yr, n_x, sum_only));
    return rcpp_result_gen;
END_RCPP
}
// em_dirichlet_wt
Eigen::VectorXd em_dirichlet_wt(const Eigen::VectorXd curr, const Eigen::VectorXi Y, const Eigen::VectorXi X, const Eigen::VectorXd wt, const Eigen::MatrixXd p_rxs, const Eigen::MatrixXd prior_yr, int n_x);
RcppExport SEXP _birdie_em_dirichlet_wt(SEXP currSEXP, SEXP YSEXP, SEXP XSEXP, SEXP wtSEXP, SEXP p_rxsSEXP, SEXP prior_yrSEXP, SEXP n_xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type curr(currSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXi >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXi >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type p_rxs(p_rxsSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type prior_yr(prior_yrSEXP);
    Rcpp::traits::input_parameter< int >::type n_x(n_xSEXP);
    rcpp_result_gen = Rcpp::wrap(em_dirichlet_wt(curr, Y, X, wt, p_rxs, prior_yr, n_x));
    return rcpp_result_gen;
END_RCPP
}
// gibbs_dir_step
Eigen::VectorXd gibbs_dir_step(const Eigen::VectorXi Y, const Eigen::VectorXi X, const Eigen::VectorXd wt, const Eigen::MatrixXd p_ryxs, const Eigen::MatrixXd prior_yr, int n_x);
RcppExport SEXP _birdie_gibbs_dir_step(SEXP YSEXP, SEXP XSEXP, SEXP wtSEXP, SEXP p_ryxsSEXP, SEXP prior_yrSEXP, SEXP n_xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXi >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXi >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type p_ryxs(p_ryxsSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type prior_yr(prior_yrSEXP);
    Rcpp::traits::input_parameter< int >::type n_x(n_xSEXP);
    rcpp_result_gen = Rcpp::wrap(gibbs_dir_step(Y, X, wt, p_ryxs, prior_yr, n_x));
    return rcpp_result_gen;
END_RCPP
}
// resid_mult
Eigen::VectorXd resid_mult(const Eigen::VectorXd m_coef, const Eigen::VectorXi idxs, const Eigen::MatrixXd r_probs, int k, int n_k);
RcppExport SEXP _birdie_resid_mult(SEXP m_coefSEXP, SEXP idxsSEXP, SEXP r_probsSEXP, SEXP kSEXP, SEXP n_kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type m_coef(m_coefSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXi >::type idxs(idxsSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type r_probs(r_probsSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type n_k(n_kSEXP);
    rcpp_result_gen = Rcpp::wrap(resid_mult(m_coef, idxs, r_probs, k, n_k));
    return rcpp_result_gen;
END_RCPP
}
// safeexpoffset
Eigen::MatrixXd safeexpoffset(const Eigen::MatrixXd Y);
RcppExport SEXP _birdie_safeexpoffset(SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(safeexpoffset(Y));
    return rcpp_result_gen;
END_RCPP
}
// gibbs_me
Eigen::MatrixXd gibbs_me(int iter, int warmup, const Eigen::VectorXi& S, const Eigen::VectorXi& GZ, const Eigen::MatrixXd& M_sr, const Eigen::MatrixXd& N_gzr, const Eigen::MatrixXd& alpha_gzr, const Eigen::MatrixXd& beta_sr, int cores, int verbosity);
RcppExport SEXP _birdie_gibbs_me(SEXP iterSEXP, SEXP warmupSEXP, SEXP SSEXP, SEXP GZSEXP, SEXP M_srSEXP, SEXP N_gzrSEXP, SEXP alpha_gzrSEXP, SEXP beta_srSEXP, SEXP coresSEXP, SEXP verbositySEXP) {
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
    Rcpp::traits::input_parameter< int >::type cores(coresSEXP);
    Rcpp::traits::input_parameter< int >::type verbosity(verbositySEXP);
    rcpp_result_gen = Rcpp::wrap(gibbs_me(iter, warmup, S, GZ, M_sr, N_gzr, alpha_gzr, beta_sr, cores, verbosity));
    return rcpp_result_gen;
END_RCPP
}
// mat_rcatp
Eigen::VectorXi mat_rcatp(Eigen::MatrixXd probs);
RcppExport SEXP _birdie_mat_rcatp(SEXP probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type probs(probsSEXP);
    rcpp_result_gen = Rcpp::wrap(mat_rcatp(probs));
    return rcpp_result_gen;
END_RCPP
}
// rdirichlet
VectorXd rdirichlet(const VectorXd alpha, const int m);
RcppExport SEXP _birdie_rdirichlet(SEXP alphaSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(rdirichlet(alpha, m));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_stan_fit4multinom_mod();

static const R_CallMethodDef CallEntries[] = {
    {"_birdie_calc_bayes_bisg", (DL_FUNC) &_birdie_calc_bayes_bisg, 5},
    {"_birdie_calc_bayes", (DL_FUNC) &_birdie_calc_bayes, 6},
    {"_birdie_dirichlet_map", (DL_FUNC) &_birdie_dirichlet_map, 5},
    {"_birdie_em_dirichlet", (DL_FUNC) &_birdie_em_dirichlet, 7},
    {"_birdie_em_dirichlet_wt", (DL_FUNC) &_birdie_em_dirichlet_wt, 7},
    {"_birdie_gibbs_dir_step", (DL_FUNC) &_birdie_gibbs_dir_step, 6},
    {"_birdie_resid_mult", (DL_FUNC) &_birdie_resid_mult, 5},
    {"_birdie_safeexpoffset", (DL_FUNC) &_birdie_safeexpoffset, 1},
    {"_birdie_gibbs_me", (DL_FUNC) &_birdie_gibbs_me, 10},
    {"_birdie_mat_rcatp", (DL_FUNC) &_birdie_mat_rcatp, 1},
    {"_birdie_rdirichlet", (DL_FUNC) &_birdie_rdirichlet, 2},
    {"_rcpp_module_boot_stan_fit4multinom_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4multinom_mod, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_birdie(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
