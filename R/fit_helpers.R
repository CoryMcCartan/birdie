
vcov_to_se <- function(vcov, map) {
    out = matrix(sqrt(diag(vcov)),
                 nrow=nrow(map), ncol=ncol(map), byrow=TRUE)
    rownames(out) = rownames(map)
    colnames(out) = colnames(map)
    out
}

to_array_yrx <- function(ests, est_dim) {
    aperm(array(ests, est_dim), c(2L, 1L, 3L))
}
to_array_xyr <- function(ests, est_dim) {
    aperm(array(ests, est_dim), c(3L, 2L, 1L))
}
to_array_xry <- function(ests, est_dim) {
    aperm(array(ests, est_dim), c(3L, 1L, 2L))
}
to_vec_xyr <- function(ests) {
    as.numeric(aperm(ests, c(3L, 2L, 1L)))
}




#' Control of BIRDiE Model Fitting
#'
#' Constructs control parameters for BIRDiE model fitting.
#' All arguments have defaults.
#'
#' @param max_iter The maximum number of EM iterations.
#' @param accel The acceleration algorithm to use in doing EM. The default
#'   `"squarem"` is good for most purposes, though `"anderson"` may be faster
#'   when there are few parameters or very tight tolerances. `"daarem"` is an
#'   excellent choice as well that works across a range of problems, though it
#'   requires installing the small `SQUAREM` package. `"none"` is not
#'   recommended unless other algorithms are running into numerical issues.
#' @param order The order to use in the acceleration algorithm. Interpretation
#'   varies by algorithm. Can range from 1 to 3 (default) for SQUAREM and from 1
#'   to the number of parameters for Anderson and DAAREM (default -1 allows the
#'   order to be determined by problem size).
#' @param anderson_restart Whether to use restarts in Anderson acceleration.
#' @param abstol The absolute tolerance used in checking convergence.
#' @param reltol The relative tolerance used in checking convergence.
#'   Ignored if `accel="squarem"`.
#'
#' @return A list containing the control parameters
#'
#' @examples
#' birdie.ctrl(max_iter=1000)
#'
#' @export
birdie.ctrl <- function(max_iter=1000, accel=c("squarem", "anderson", "daarem", "none"),
                        order=switch(match.arg(accel), none=0L, anderson=-1L, daarem=-1L, squarem=3L),
                        anderson_restart=TRUE,
                        abstol=1e-6, reltol=1e-6) {
    stopifnot(max_iter >= 1)
    stopifnot(abstol >= 0)
    stopifnot(reltol >= 0)

    accel = match.arg(accel)
    fn_accel = switch(accel,
                      none = accel_none,
                      anderson = accel_anderson,
                      daarem = accel_daarem,
                      squarem = accel_squarem)
    if (accel == "squarem") order = min(order, 3)

    list(max_iter = as.integer(max_iter),
         accel = fn_accel,
         order = order,
         restart = as.logical(anderson_restart),
         abstol = abstol,
         reltol = reltol)
}

# Check (and possibly create default) prior
check_make_prior <- function(prior, method, n_y, n_r) {
    if (is.null(prior)) {
        prior = matrix(1, nrow=n_y, ncol=n_r)
        if (method == "pool") {
            cli_inform("Using uniform prior for Pr(X | R)",
                       call=parent.frame())
        } else if (method == "fixef") {
            cli_inform("Using c(1+\u03B5, 1+\u03B5, ..., 1+\u03B5) prior for Pr(X | R)",
                       call=parent.frame())
            prior = matrix(1 + 100*.Machine$double.eps, nrow=n_y, ncol=n_r)
        }
    }
    if (nrow(prior) != n_y) {
        cli_abort("{.arg prior} must have the same number of rows
                  as there are levels of X", call=parent.frame())
    }
    if (ncol(prior) != n_r) {
        cli_abort("{.arg prior} must have the same number of columns
                  as there are racial groups", call=parent.frame())
    }

    prior
}

# Check predictors against theory
check_covars <- function(r_probs, covars, method) {
    if (inherits(r_probs, "bisg")) {
        if (attr(r_probs, "S_name") %in% covars) {
            cli_warn("Last name vector {.arg {attr(r_probs, 'S_name')}}
                     should not be used in BIRDiE model",
                     call=parent.frame())
        }
        if (method != "pool" &&
            length(setdiff( attr(r_probs, "GX_names"), covars)) > 0) {
            cli_warn(c("Some variables used to create BISG probabilities are
                       not used in BIRDiE model.",
                       "x"="Statistically valid inference is not guaranteed."),
                     call=parent.frame())
        }
    }
}

# Check for full interaction structure
check_full_int <- function(tt, covars) {
    int_ord = attr(tt, "order")
    int_ord = table(factor(int_ord, levels=seq_len(max(int_ord))))
    full_int = all(int_ord == choose(length(int_ord), seq_along(int_ord))) # pascal's triangle
    if (!full_int) {
        x_vars = attr(tt, "term.labels")[attr(tt, "order") == 1]
        cli_warn(c("Fixed effects (no-pooling) model is specified without full
                   interaction structure but will be fit with full interactions.",
                   "i"="For full interaction structure, use
                         `{covars[1]} ~ {paste0(x_vars, collapse=' * ')}`."),
                 call=parent.frame())
    }
}
