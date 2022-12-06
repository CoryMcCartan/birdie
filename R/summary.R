#' Compute Racial Disparities from Model Estimates
#'
#' This function lets you easily compute differences in conditional expectations
#' between all pairs of specified racial groups.
#'
#' @param x A `birdie` model object.
#' @param subgroup If `TRUE`, return subgroup-level (rather than marginal)
#'   disparity estimates.
#' @param races A character vector of racial groups to compute disparities for.
#'   The special value `TRUE`, the default, computes disparities for all racial
#'   groups.
#'
#' @return A tibble containing a row with every possible disparity for the
#'   specified `races`, which are identified by columns `race_1` and `race_2`.
#'   The reported disparity is `estimate_1 - estimate_2`.
#'
#' @examples
#' data(pseudo_vf)
#' r_probs = bisg(~ nm(last_name) + zip(zip), data=pseudo_vf)
#' fit = birdie(r_probs, turnout ~ 1, data=pseudo_vf)
#'
#' disparities(fit)
#' disparities(fit, races=c("white", "black"))
#'
#' @concept estimators
#' @export
disparities <- function(x, subgroup=FALSE, races=TRUE) {
    if (isTRUE(races)) {
        races = colnames(x$map)
    }

    td = filter(tidy(x, subgroup=subgroup), .data$race %in% races)
    gx_names = setdiff(colnames(td), c("race", "estimate"))

    combos = expand.grid(race_1=races, race_2=races) %>%
        filter(.data$race_1 != .data$race_2) %>%
        as_tibble()

    left_join(combos, td, by=c("race_1"="race")) %>%
        left_join(td, by=c("race_2"="race", gx_names),
                  suffix=c("_1", "_2")) %>%
        mutate(disparity = .data$estimate_1 - .data$estimate_2)
}
