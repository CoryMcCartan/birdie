
#' Estimate the joint distribution of R and X
#'
#' @param r_probs a matrix data frame of race probabilities
#' @param x the variable X to tabulate against
#' @param method the estimation method: a weighted-mean estimator, a
#'   thresholding estimator, or a multiple-imputation estimator.
#' @param prefix how to select the race probability columns from `r_probs`, if
#'   if is a data frame
#' @param n_mi if `method == "mi"`, how many imputations to produce.
#'
#' @returns a matrix
#' @export
calc_joint_bisgz = function(r_probs, x, method=c("weight", "thresh", "mi", "ols"),
                            prefix="pr_", n_mi=8) {
    x = as.factor(x)
    N = length(x)
    if (!is.matrix(r_probs)) {
        r_probs = as.matrix(dplyr::select(r_probs, starts_with(prefix)))
        colnames(r_probs) = substring(colnames(r_probs), nchar(prefix)+1L)
    }
    stopifnot(nrow(r_probs) == N)

    method = match.arg(method)
    if (method == "weight") {
        out = lapply(levels(x), function(l) colMeans(r_probs * (x == l)))
        out = do.call(rbind, out)
        rownames(out) = levels(x)
    } else if (method == "thresh") {
        r_est = factor(max.col(r_probs), levels=seq_len(ncol(r_probs)))
        out = table(x, r_est) / N
        names(dimnames(out)) = NULL
        colnames(out) = colnames(r_probs)
    } else if (method == "mi") {
        n_race = ncol(r_probs)
        out = matrix(0, nrow=nlevels(x), ncol=n_race)
        for (i in seq_len(n_mi)) {
            r_est = apply(r_probs, 1, function(x) sample(n_race, 1, prob=x))
            out = out + table(x, r_est) / N / n_mi
        }
        names(dimnames(out)) = NULL
        colnames(out) = colnames(r_probs)
    } else if (method == "ols") {
        out = lapply(levels(x), function(l) {
            lm.fit(r_probs, x == l)$coefficients
        })
        out = do.call(rbind, out)
        out = out %*% diag(colMeans(r_probs))
        rownames(out) = levels(x)
    }
    out
}

#' Calculate a posterior quantile of the joint distribution of R and X
#'
#' @param fit a `fit_raceproxy` object (the output of `model_race`)
#' @param which if `condition` was used in `model_race`, which estimates to extract.
#' @param q the quantile
#' @param p_r a vector containing the marginal probabilities for each value of
#'   `R`. Defaults to the demographics of the US.
#'
#' @returns a matrix
#' @export
calc_joint_model = function(fit, which="global", q=0.5,
                            p_r=c(white=0.630, black=0.121, hisp=0.173,
                                  asian=0.0478, aian=0.0072, other=0.0210)) {
    out = apply(fit$draws[[which]], 2:3, function(x) quantile(x, q)) %*% diag(p_r)
    colnames(out) = fit$r_lev
    rownames(out) = fit$x_lev
    out
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
eval_joints = function(tgt, metric=c("tv", "tv_col", "tv_row", "mad", "rmse"), ...) {
    score_fn = function(x) cli_abort("Metric {.val {metric}} not recognized.")
    metric = match.arg(metric)

    n_out = 1
    if (metric == "tv") {
        score_fn = function(x) sum(abs(tgt - x)) / 2
    } else if (metric == "tv_col") {
        score_fn = function(x) list(colSums(abs(tgt - x)) / 2)
        n_out = ncol(tgt)
    } else if (metric == "tv_row") {
        score_fn = function(x) list(rowSums(abs(tgt - x)) / 2)
        n_out = nrow(tgt)
    } else if (metric == "mad") {
        score_fn = function(x) mean(abs(tgt - x))
    } else if (metric == "rmse") {
        score_fn = function(x) sqrt(mean((tgt - x)^2))
    }

    joints = rlang::list2(...) %>%
        #vapply(score_fn, numeric(n_out)) %>%
        sapply(score_fn)# %>%
    if (n_out == 1) {
        joints = sort(joints)
    }
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
