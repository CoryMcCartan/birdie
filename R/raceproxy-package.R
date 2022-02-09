#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
#' @importFrom cli cli_inform cli_warn cli_abort
#' @import dplyr
#' @importFrom reticulate py
## usethis namespace: end
NULL

.onLoad <- function(libname, pkgname) {
    reticulate::configure_environment(pkgname)
}
