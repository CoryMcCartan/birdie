#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
#' @importFrom cli cli_inform cli_warn cli_abort
#' @import dplyr
#' @importFrom Rcpp sourceCpp
#' @importFrom rstan sampling
#' @importFrom reticulate py
#' @useDynLib raceproxy, .registration = TRUE
## usethis namespace: end
NULL

.onLoad <- function(libname, pkgname) {
    reticulate::configure_environment(pkgname)
}
