#' Sensitivity Analysis for BIRDiE Models
#'
#' @param x A `birdie` model object.
#' @param r_probs A data frame of BISG probabilities, from [bisg()].
#' @param contrast A named numeric vector specifying the estimate, disparity, or
#'   other contrast to compute a sensitivity bound for. Must have length no more
#'   than the number of racial groups. Missing racial groups will have a
#'   coefficient of zero. For example, to get a bound on disparities between
#'   White and Black individuals, use `contrast=c(white=1, black=-1, hisp=0, ...)`.
#'
#' @return A vector containing the sensitivity bound factor for the selected
#'   contrast, for each level of the outcome variable. This factor should be
#'   multiplied by the assumed partial \eqn{R^2} of the outcome residuals on
#'   surnames in order to produce a bound on the estimate of the bias
#'
#' @export
sens_bound <- function(x, r_probs, contrast) {
    if (!inherits(x, "birdie"))
        cli_abort("{.arg x} must be a {.cls birdie} object.")
    if (!inherits(r_probs, "bisg"))
        cli_abort("{.arg r_probs} must be a {.cls bisg} object from {.fn bisg}.")

    races = colnames(x$map)
    contr = rep_along(races, 0)
    names(contr) = races
    contr[names(contrast)] = contrast

    if (is.null(attr(r_probs, "p_rgx"))) {
        cli_abort(c("Missing Pr(R | G, X) table; cannot generate bounds.",
                    ">"="Generate your {.fn bisg} predictions with
                        {.arg save_rgx=TRUE}."))
    }
    eps = 5e-4
    p_rgxs = as.matrix(r_probs) * (1 - eps) + eps
    p_rgx = (attr(r_probs, "p_rgx") * (1 - eps) + eps)[attr(r_probs, "gx"), ]

    e_alph2 = sqrt(sum(contr^2 * colMeans(1/p_rgxs - 1/p_rgx)))
    e_res = colMeans(residuals.birdie(x, x_only=TRUE)^2)

    sqrt(e_res * e_alph2)
}
