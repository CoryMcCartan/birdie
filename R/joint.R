
#' Evaluate the quality of a joint distribution approximation
eval_joints = function(tgt, metric=c("tv", "mad", "rmse"), ...) {
    score_fn = function(x) cli_abort("Metric {.val {metric}} not recognized.")
    if (metric == "tv") {
        score_fn = \(x) sum(abs(tgt - x)) / 2
    } else if (metric == "mad") {
        score_fn = \(x) mean(abs(tgt - x))
    } else if (metric == "rmse") {
        score_fn = \(x) sqrt(mean((tgt - x)^2))
    }

    joints = rlang::list2(...) %>%
        vapply(score_fn, numeric(1)) %>%
        sort()
    tibble(method = names(joints),
           "{toupper(metric)}" := joints)
}

#' Rake a joint distribution so that the margins match provided vectors
#'
#' @param est the starting joint distribution
#' @param row the desired row sums / margins
#' @param col the desired column sums / margins
#'
#' @returns a new matrix `est`
#'
#' @export
rake = function(est, row=rowSums(est), col=colSums(est)) {
    for (i in 1:5) {
        est = est %*% diag(col / colSums(est))
        est = diag(row / rowSums(est)) %*% est
    }
    est
}
