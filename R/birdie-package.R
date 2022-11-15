#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
#' @import dplyr
#' @importFrom rlang expr enquo enquos sym syms .data := as_name as_label eval_tidy
#' @importFrom tidyselect eval_select
#' @importFrom cli cli_inform cli_warn cli_abort
#' @importFrom Rcpp sourceCpp
#' @importFrom reticulate py
#' @importFrom stats na.omit quantile lm.fit
#' @importFrom utils tail
#' @useDynLib birdie, .registration = TRUE
## usethis namespace: end
NULL

# module with python code
py_code = NULL

.onLoad <- function(libname, pkgname) {
    reticulate::configure_environment(pkgname)
    py_path = system.file("py", package="birdie")
    py_code <<- reticulate::import_from_path("birdie", path=py_path, delay_load=FALSE)
}