#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
#' @import dplyr
#' @importFrom rlang expr enquo enquos sym syms .data :=
#' @importFrom rlang as_name as_label eval_tidy rep_along f_lhs f_rhs
#' @importFrom stringr str_c str_detect str_starts str_remove_all str_replace_all str_remove
#' @importFrom cli cli_inform cli_warn cli_abort
#' @importFrom Rcpp sourceCpp
#' @importFrom stats na.omit na.fail .lm.fit fitted simulate median cov
#' @importFrom stats terms get_all_vars update.formula as.formula model.frame model.matrix
#' @importFrom utils tail unstack relist head
#' @import methods
#' @rawNamespace useDynLib(birdie, .registration=TRUE)
## usethis namespace: end
NULL

Rcpp::loadModule("stan_fit4multinom_mod", what=TRUE)
