// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "raceproxy_types.h"
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// gibbs_race
mat gibbs_race(int iter, int warmup, const uvec& X, const uvec& S, const uvec& GZ, const mat& M_sr, const mat& N_gzr, const mat& alpha_sr, const mat& beta_gzr, int verbosity);
RcppExport SEXP _raceproxy_gibbs_race(SEXP iterSEXP, SEXP warmupSEXP, SEXP XSEXP, SEXP SSEXP, SEXP GZSEXP, SEXP M_srSEXP, SEXP N_gzrSEXP, SEXP alpha_srSEXP, SEXP beta_gzrSEXP, SEXP verbositySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< int >::type warmup(warmupSEXP);
    Rcpp::traits::input_parameter< const uvec& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const uvec& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const uvec& >::type GZ(GZSEXP);
    Rcpp::traits::input_parameter< const mat& >::type M_sr(M_srSEXP);
    Rcpp::traits::input_parameter< const mat& >::type N_gzr(N_gzrSEXP);
    Rcpp::traits::input_parameter< const mat& >::type alpha_sr(alpha_srSEXP);
    Rcpp::traits::input_parameter< const mat& >::type beta_gzr(beta_gzrSEXP);
    Rcpp::traits::input_parameter< int >::type verbosity(verbositySEXP);
    rcpp_result_gen = Rcpp::wrap(gibbs_race(iter, warmup, X, S, GZ, M_sr, N_gzr, alpha_sr, beta_gzr, verbosity));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_stan_fit4simple_additive_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4simple_nonparam_mod();

static const R_CallMethodDef CallEntries[] = {
    {"_raceproxy_gibbs_race", (DL_FUNC) &_raceproxy_gibbs_race, 10},
    {"_rcpp_module_boot_stan_fit4simple_additive_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4simple_additive_mod, 0},
    {"_rcpp_module_boot_stan_fit4simple_nonparam_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4simple_nonparam_mod, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_raceproxy(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
