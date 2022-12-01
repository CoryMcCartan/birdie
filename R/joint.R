
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
        cli_inform("Use {.fn calc_joint_bisgz_ols} for the unbiased
                    poststratified OLS estimator. This function implements
                    a simpler OLS estimator that may be biased.",
                   .frequency="once", .frequency_id="calc_joint")
        out = lapply(levels(x), function(l) {
            lm.fit(r_probs, x == l)$coefficients
        })
        out = do.call(rbind, out)
        out = out %*% diag(colMeans(r_probs))
        rownames(out) = levels(x)
        colnames(out) = colnames(r_probs)
    }
    out
}

#' Estimate the joint distribution of R and X
#'
#' @param r_probs a matrix data frame of race probabilities
#' @param x the variable X to tabulate against
#' @param gz values of GZ variables. Usually exported as part of
#'   [predict_race_sgz()]
#' @param p_gzr a table giving Pr(GZ|R) for use in post-stratification. Columns
#'   should sum to 1.
#' @param prefix how to select the race probability columns from `r_probs`, if
#'   if is a data frame
#' @param truncate if `TRUE`, truncate results to lie in \[0, 1\].
#'
#' @returns a matrix
#' @export
calc_joint_bisgz_ols = function(r_probs, x, gz=attr(r_probs, "gz"),
                                p_gzr=attr(r_probs, "p_gzr"),
                                prefix="pr_", truncate=TRUE) {
    x = as.factor(x)
    N = length(x)
    gz = as.factor(gz)
    if (!is.matrix(r_probs)) {
        r_probs = as.matrix(dplyr::select(r_probs, starts_with(prefix)))
        colnames(r_probs) = substring(colnames(r_probs), nchar(prefix)+1L)
    }
    stopifnot(nrow(r_probs) == N)

    yy = do.call(cbind, lapply(levels(x), function(xval) {
        as.numeric(x == xval)
    }))
    ests_local = do.call(rbind, lapply(levels(gz), function(lvl) {
        idx = which(gz == lvl)
        tryCatch(
            as.numeric(.lm.fit(r_probs[idx, ], yy[idx, ])$coefficients),
        error = function(e) rep(0, 6))
    }))
    if (isTRUE(truncate)) {
        ests_local[ests_local < 0] = 0
        ests_local[ests_local > 1] = 1
    }
    out = matrix(nrow=nlevels(x), ncol=ncol(r_probs))
    for (i in seq_len(nlevels(x))) {
        out[i, ] = colSums(p_gzr * ests_local[, 1:6 + 6*(i-1)], na.rm=TRUE)
    }


    out = out %*% diag(colMeans(r_probs))
    rownames(out) = levels(x)
    colnames(out) = colnames(r_probs)

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
        # score_fn = function(x) sum(abs(tgt - x)) / 2
        score_fn = function(x) sum(abs(tgt[, 1:5] - x[, 1:5])) / 2
    } else if (metric == "tv_col") {
        score_fn = function(x) list(colSums(abs(tgt - x)) / 2 / colSums(tgt))
        n_out = ncol(tgt)
    } else if (metric == "tv_row") {
        score_fn = function(x) list(rowSums(abs(tgt - x)) / 2 / rowSums(tgt))
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
#' @param iter the number of raking iterations
#'
#' @returns a new matrix `est`
#'
#' @export
rake = function(est, row=rowSums(est), col=colSums(est), iter=5) {
    for (i in seq_len(iter)) {
        est = est %*% diag(col / colSums(est))
        est = diag(row / rowSums(est)) %*% est
    }
    est
}
