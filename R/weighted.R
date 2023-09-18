#' Calculate Weighted Estimate of Outcomes By Race
#'
#' Calculates the "standard" weighted estimator of conditional distributions of
#' an outcome variable \eqn{Y} by race \eqn{R}, using BISG probabilities.  This
#' estimator, while commonly used, is only appropriate if \eqn{Y \perp R \mid X, S},
#' where \eqn{S} and \eqn{X} are the last names and covariates (possibly
#' including geography) used in making the BISG probabilities. In most cases
#' this assumption is not plausible and [birdie()] should be used instead. See
#' the references below for more discussion as to selecting the right estimator.
#'
#' @param r_probs A data frame or matrix of BISG probabilities, with one row per
#'   individual. The output of [bisg()] can be used directly here.
#' @param formula A two-sided formula object describing the estimator structure.
#'   The left-hand side is the outcome variable, which must be discrete.
#'   Subgroups for which to calculate estimates may be specified by adding
#'   covariates on the right-hand side. Subgroup estimates are available with
#'   `coef(..., subgroup=TRUE)` and `tidy(..., subgroup=TRUE)`.
#' @param data An optional data frame containing the variables named in `formula`.
#' @param prefix If `r_probs` is a data frame, the columns containing racial
#'   probabilities will be selected as those with names starting with `prefix`.
#'   The default will work with the output of [bisg()].
#' @param se_boot The number of bootstrap replicates to use to compute
#'   approximate standard errors for the estimator. When there are fewer than
#'   1,000 individuals or 100 or fewer replicates, a Bayesian bootstrap is used
#'   instead (i.e., weights are drawn from a \eqn{\text{Dirichlet}(1, 1, ...,
#'   1)} distribution, which produces more reliable estimates.
#'
#' @return An object of class `est_weighted`, inheriting from
#'   [`birdie`][birdie::birdie-class], for which many methods are available. The
#'   model estimates may be accessed with `coef()`.
#'
#' @references
#' McCartan, C., Fisher, R., Goldin, J., Ho, D., & Imai, K. (2022).
#' Estimating Racial Disparities when Race is Not Observed.
#' Available at \url{https://arxiv.org/abs/2303.02580}.
#'
#' @examples
#' data(pseudo_vf)
#'
#' r_probs = bisg(~ nm(last_name) + zip(zip), data=pseudo_vf)
#'
#' # Process zip codes to remove missing values
#' pseudo_vf$zip = proc_zip(pseudo_vf$zip)
#'
#' est_weighted(r_probs, turnout ~ 1, data=pseudo_vf)
#'
#' est = est_weighted(r_probs, turnout ~ zip, data=pseudo_vf)
#' tidy(est, subgroup=TRUE)
#'
#' @concept estimators
#' @export
est_weighted <- function(r_probs, formula, data=NULL, prefix="pr_", se_boot=0) {
    Y_vec = eval_tidy(f_lhs(formula), data)
    tt = terms(formula, keep.order=TRUE)
    covars = all.vars(tt)

    # set up race probability matrix
    if (is.matrix(r_probs)) {
        p_rxs = r_probs
    } else if (is.data.frame(r_probs)) {
        p_rxs = as.matrix(select(r_probs, starts_with(prefix)))
    } else {
        cli_abort("{.arg r_probs} must be a matrix or data frame.")
    }

    # check types
    if (!check_discrete(Y_vec))
        cli_abort("Response variable must be a character or factor with no missing values.")
    if (!is.null(data) && nrow(data) != nrow(r_probs))
        cli_abort("{.arg data} and {.arg r_probs} must have the same number of rows.")

    Y_vec = as.factor(Y_vec)
    n_y = nlevels(Y_vec)
    n_r = ncol(p_rxs)
    d_model = model.frame(formula, data=data, na.action=na.fail)[-1]
    Y = as.integer(Y_vec)

    # create unique group IDs
    X = to_unique_ids(d_model)
    idx_sub = vctrs::vec_unique_loc(X)
    tbl_gx = d_model[idx_sub, , drop=FALSE]
    n_x = max(X)
    est_dim = c(n_r, n_y, n_x)
    ones_mat = matrix(1, nrow=n_y, ncol=n_r)

    # compute estimates
    est_glb = dirichlet_map(Y, rep_along(Y, 1), p_rxs, ones_mat, 1) %>%
        matrix(n_y, n_r, byrow=TRUE)
    est_sub = dirichlet_map(Y, X, p_rxs, ones_mat, n_x) %>%
        to_array_yrx(est_dim)

    # add names
    colnames(est_glb) = stringr::str_sub(colnames(p_rxs), nchar(prefix)+1L)
    rownames(est_glb) = levels(Y_vec)
    dimnames(est_sub) = c(dimnames(est_glb), list(tbl_gx_names(tbl_gx)))

    # format p_ryxs
    p_ryxs = as_tibble(p_rxs)
    if (inherits(r_probs, "bisg")) {
        p_ryxs = reconstruct.bisg(p_ryxs, r_probs)
    }

    # bootstrap
    if (se_boot > 0) {
        boot_ests = boot_wtd(se_boot, Y, n_y, p_rxs)
        vcov = cov(t(boot_ests))
    }

    # output
    attr(tt, ".Environment") = NULL # save space

    structure(list(
        map = est_glb,
        map_sub = est_sub,
        p_ryxs = p_ryxs,
        vcov = if (se_boot > 0) vcov else NULL,
        se = if (se_boot > 0) vcov_to_se(vcov, est_glb) else NULL,
        N = length(Y_vec),
        tbl_gx = as_tibble(tbl_gx),
        vec_gx = X,
        y_name = covars[1],
        y = Y_vec,
        prefix = prefix,
        entropy = list(pre = median(entropy(p_rxs)),
                       post = median(entropy(p_ryxs))),
        call = match.call()
    ), class=c("est_weighted", "birdie"))
}


boot_wtd <- function(R=10, Y, n_y, p_rxs) {
    N = length(Y)
    n_r = ncol(p_rxs)
    ones = rep_along(Y, 1)
    ones_mat = matrix(1, nrow=n_y, ncol=n_r)

    out = matrix(nrow=n_r*n_y, ncol=R)

    if (N > 1000 && R > 100) {
        mk_wt = function() tabulate(sample.int(N, N, replace=TRUE), N)
    } else { # more computationally intensive but smoother
        mk_wt = function() N * as.numeric(rdirichlet(1, ones))
    }

    cli::cli_progress_bar("Bootstrapping", total=R)
    for (i in seq_len(R)) {
        wt = mk_wt()

        out[, i] = dirichlet_map(Y, ones, p_rxs*wt, ones_mat, 1)

        cli::cli_progress_update()
    }
    cli::cli_progress_done()

    out
}


#' @describeIn est_weighted Print a summary of the model fit.
#' @param object,x An object of class `est_weighted`.
#' @param ... Additional arguments to generic methods (ignored).
#' @export
print.est_weighted <- function(x, ...) {
    cli::cli_text("Weighted estimator")
    cli::cat_line("Formula: ", deparse(x$call$formula))
    cli::cat_line("   Data: ", deparse(x$call$data))
    cli::cli_text("Number of obs: {comma(x$N)};
                  groups: {comma(dim(x$map_sub)[3])}")

    cli::cli_text("Estimated distribution:")
    m = round(x$map, 3)
    print(m)
}

#' @describeIn est_weighted Print a more detailed summary of the model fit.
#' @export
summary.est_weighted <- function(object, ...) {
    cli::cli_text("Weighted estimator")
    cli::cat_line("Formula: ", deparse(object$call$formula))
    cli::cat_line("   Data: ", deparse(object$call$data))
    cat("\n")

    cli::cli_text("Number of observations: {comma(object$N)}")
    cli::cli_text("Number of groups: {comma(dim(object$map_sub)[3])}")
    cat("\n")

    p_r = colMeans(object$p_ryxs)
    cli::cat_line("BISG entropy decrease from marginal race distribution:")
    print(summary(entropy(p_r) - entropy(object$p_ryxs)))
    cat("\n")

    cli::cli_text("Estimated outcome-by-race distribution:")
    m = round(object$map, 3)
    print(m)

    invisible(object)
}

