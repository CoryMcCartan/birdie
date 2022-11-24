
#' @export
print.birdie <- function(x, ...) {
    cli::cli_text("A {.pkg BIRDiE} model fit with
                  {format(x$N, big.mark=',')} observations")
    cat("\n")

    cli::cli_text("Estimated distribution:")
    cat("\n")
    m = round(x$map, 3)
    colnames(m) = x$r_lev
    rownames(m) = x$x_lev
    print(m)
}

