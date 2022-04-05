#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
#' @import dplyr
#' @importFrom cli cli_inform cli_warn cli_abort
#' @importFrom Rcpp sourceCpp
#' @importFrom reticulate py
#' @importFrom stats na.omit quantile
#' @useDynLib raceproxy, .registration = TRUE
## usethis namespace: end
NULL

# module with python code
py_code = NULL

.onLoad <- function(libname, pkgname) {
    reticulate::configure_environment(pkgname)
    py_path = system.file("py", package="raceproxy")
    py_code <<- reticulate::import_from_path("raceproxy", path=py_path, delay_load=FALSE)
}
