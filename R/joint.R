
#' Calculate the joint distribution of R and X
#'
#' @param r_probs a matrix data frame of race probabilities
#' @param x the variable X to tabulate against
#' @param prefix how to select the race probability columns from `r_probs`, if
#'   if is a data frame
#'
#' @returns a matrix
#' @export
calc_joint_bisgz = function(r_probs, x, prefix="pr_") {
    x = as.factor(x)
    if (!is.matrix(r_probs)) {
        r_probs = as.matrix(dplyr::select(r_probs, starts_with(prefix)))
        colnames(r_probs) = substring(colnames(r_probs), nchar(prefix)+1L)
    }

    lapply(levels(x), function(.) colMeans(r_probs * (x == .))) %>%
        do.call(rbind, .) %>%
        `rownames<-`(levels(x))
}

#' Calculate a posterior quantile of the joint distribution of R and X
#'
#' @param draws a 3-dim array with the first dimension representing draws
#' @param q the quantile
#' @param p_r a vector containing the marginal probabilities for each value of
#'   `R`. Defaults to the demographics of the US.
#'
#' @returns a matrix
#' @export
calc_joint_model = function(draws, q=0.5,
                            p_r=c(white=0.615, black=0.123, hisp=0.176, asian=0.053, other=0.034)) {
    apply(draws, 2:3, \(x) quantile(x, q)) %*% diag(p_r) %>%
        `colnames<-`(names(p_r))
}


#' Evaluate the quality of a joint distribution approximation
#'
#' @param tgt the true joint distribution
#' @param metric how to evaluate the joint distributions: total variation
#'   (`tv`), mean absolute deviation (`mad`), or root mean-square error
#'   (`rmse`).
#' @param ... named joint distributions.
#'
#' @return a tibble with a row for every argument in `...`
#' @export
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
