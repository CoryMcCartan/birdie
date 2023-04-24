
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
#'   requires installing the small `daarem` package. `"none"` is not
#'   recommended unless other algorithms are running into numerical issues.
#'   See the references below for details on these schemes.
#' @param order The order to use in the acceleration algorithm. Interpretation
#'   varies by algorithm. Can range from 1 (default) to 3 for SQUAREM and from 1
#'   to the number of parameters for Anderson and DAAREM (default -1 allows the
#'   order to be determined by problem size).
#' @param anderson_restart Whether to use restarts in Anderson acceleration.
#' @param abstol The absolute tolerance used in checking convergence.
#' @param reltol The relative tolerance used in checking convergence.
#'   Ignored if `accel = "squarem"` or `"daarem"`.
#'
#' @return A list containing the control parameters.
#'
#' @references
#' Varadhan, R., & Roland, C. (2004). Squared extrapolation methods (SQUAREM): A
#' new class of simple and efficient numerical schemes for accelerating the
#' convergence of the EM algorithm.
#'
#' Walker, H. F., & Ni, P. (2011). Anderson acceleration for fixed-point
#' iterations. SIAM Journal on Numerical Analysis, 49(4), 1715-1735.
#'
#' Henderson, N. C., & Varadhan, R. (2019). Damped Anderson acceleration with
#' restarts and monotonicity control for accelerating EM and EM-like algorithms.
#' Journal of Computational and Graphical Statistics, 28(4), 834-846.
#'
#' @examples
#' str(birdie.ctrl(max_iter=100))
#'
#' @concept estimators
#' @export
birdie.ctrl <- function(max_iter=1000, accel=c("squarem", "anderson", "daarem", "none"),
                        order=switch(match.arg(accel), none=0L, anderson=-1L, daarem=-1L, squarem=1L),
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

vcov_to_se <- function(vcov, map) {
    out = matrix(sqrt(diag(vcov)),
                 nrow=nrow(map), ncol=ncol(map), byrow=TRUE)
    rownames(out) = rownames(map)
    colnames(out) = colnames(map)
    out
}

tbl_gx_names <- function(tbl_gx) {
    if (nrow(tbl_gx) <= 1) {
        NULL
    } else if (ncol(tbl_gx) == 1) {
        tbl_gx[[1]]
    } else {
        for (i in seq_len(ncol(tbl_gx))) {
            if (is.numeric(tbl_gx[[i]])) {
                tbl_gx[[i]] = abbreviate(tbl_gx[[i]], 1)
            } else {
                tbl_gx[[i]] = as.character(tbl_gx[[i]])
            }
        }
        Reduce(function(x, y) str_c(x, "/", y), tbl_gx)
    }
}

to_array_yrx <- function(ests, est_dim) {
    aperm(array(ests, est_dim), c(2L, 1L, 3L))
}
to_array_xyr <- function(ests, est_dim) {
    aperm(array(ests, est_dim), c(3L, 2L, 1L))
}
to_ests_vec <- function(par_l, n_y, n_r, n_x) {
    out = array(dim=c(n_r, n_x, n_y))
    for (j in seq_along(par_l)) {
        out[j, , ] = exp(par_l[[j]]$lsft)
    }
    as.numeric(aperm(out, c(1L, 3L, 2L)))
}


# Find "|" in formula
logi_ranef <- function(formula) {
    if (!inherits(formula, "terms")) {
        formula = terms(formula)
    }
    term_labs = attr(formula, "term.labels")
    str_detect(term_labs, stringr::fixed("|"))
}
count_ranef <- function(formula) {
    sum(logi_ranef(formula))
}
re_terms <- function(formula) {
    if (!inherits(formula, "terms")) {
        formula = terms(formula)
    }
    term_labs = attr(formula, "term.labels")
    term_labs[logi_ranef(formula)]
}
remove_ranef <- function(formula) {
    if (!inherits(formula, "terms")) {
        formula = terms(formula)
    }
    re = re_terms(formula)
    if (length(re) == 0) return(formula)
    re_adj = str_c("(", re, ")", collapse=" - ")
    update.formula(formula, str_c("~ . -", re_adj))
}



# Check (and possibly create default) prior
check_make_prior <- function(prior, model, n_y, n_r) {
    if (is.null(prior)) {
        if (model == "dir") {
            cli_inform("Using c(1+\u03B5, 1+\u03B5, ..., 1+\u03B5) prior for Pr(X | R)",
                       .frequency="once", .frequency_id="birdie_prior_dir",
                       call=parent.frame())
            prior = list(
                alpha = matrix(1 + 100*.Machine$double.eps, nrow=n_y, ncol=n_r)
            )
        } else if (model == "mmm") {
            prior = list(
                scale_sigma = 0.2,
                scale_beta = 1.0
            )

            cli_inform(c("Using default prior for Pr(X | R):",
                         ">"="Prior scale on fixed effects coefficients:
                              {format(prior$scale_beta, nsmall=1)}",
                         ">"="Prior mean of random effects standard deviation:
                              {format(prior$scale_sigma, nsmall=2)}"),
                       .frequency="once", .frequency_id="birdie_prior_mmm",
                       call=parent.frame())
        }
    }

    if (model == "dir") {
        if (!"alpha" %in% names(prior) ||
                !is.numeric(prior$alpha) || !is.matrix(prior$alpha) ||
                any(is.na(prior$alpha))) {
            cli_abort(c("With {.arg model=\"dir\"}, {.arg prior} must have an entry
                        {.code alpha} which is a numeric matrix.",
                        "i"="See {.fn birdie::birdie} for details."),
                        call=parent.frame())
        }
        if (nrow(prior$alpha) != n_y) {
            cli_abort("{.arg prior$alpha} must have the same number of rows
                      as there are levels of X", call=parent.frame())
        }
        if (ncol(prior$alpha) != n_r) {
            cli_abort("{.arg prior$alpha} must have the same number of columns
                      as there are racial groups", call=parent.frame())
        }
        if (any(prior$alpha < 0)) {
            cli_abort("{.arg prior$alpha} must have nonnegative entries", call=parent.frame())
        }
        if (any(prior$alpha <= 1)) {
            cli_warn("A {.arg prior$alpha} with entries that are not
                     strictly greater than 1 may lead to numerical
                     issues.", call=parent.frame())
        }
    } else if (model == "mmm") {
        if (!all(c("scale_sigma", "scale_beta") %in% names(prior)) ||
                !is.numeric(prior$scale_beta) || length(prior$scale_beta) != 1 ||
                !is.numeric(prior$scale_sigma) || length(prior$scale_sigma) != 1) {
            cli_abort(c("With {.arg model=\"mmm\"}, {.arg prior} must have two
                        scalar entries {.code scale_beta} and {.code scale_sigma}.",
                        "i"="See {.fn birdie::birdie} for details."),
                      call=parent.frame())
        }
    }

    prior
}

# Check predictors against theory
check_covars <- function(r_probs, covars, model) {
    if (inherits(r_probs, "bisg")) {
        if (attr(r_probs, "S_name") %in% covars) {
            cli_warn("Last name vector {.arg {attr(r_probs, 'S_name')}}
                     should not be used in BIRDiE model.",
                     call=parent.frame())
        }
        unused = setdiff(attr(r_probs, "GX_names"), covars)
        if (length(covars) > 1 &&
            length(unused) > 0) {
            cli_warn(c("Some variables used to create BISG probabilities are
                       not used in BIRDiE model.",
                       "Missing: {.code {unused}}.",
                       "x"="Statistically valid inference is not guaranteed."),
                     call=parent.frame())
        }
    }
}

check_model <- function(model, tt, covars, full_int, se_boot) {
    if (model == "dir") {
        if (!full_int) {
            x_vars = attr(tt, "term.labels")[attr(tt, "order") == 1]
            cli_abort(c("Fixed effects (no-pooling) model is specified without full
                        interaction structure.",
                        "i"="Models without full interaction structure must be fit
                             using {.arg model = \"mmm\"}.",
                        ">"="For full interaction structure, use
                              `{covars[1]} ~ {paste0(x_vars, collapse=' * ')}`."),
                      call=parent.frame())
        }
        if (count_ranef(tt) > 0) {
            cli_abort("Models with random effects must be fit
                             using {.arg model = \"mmm\"}.",
                      call=parent.frame())
        }
    } else if (model == "mmm") {
        re = re_terms(tt)

        if (length(re) > 1) {
            cli_abort(c("{.pkg birdie} currently supports at most one random intercept term.",
                        "x"="Found: {.code {re}}."),
                      call=parent.frame())
        }

        if (length(re) == 1) {
            if (!str_starts(re, "1 | ") || str_detect(re, "[*:+]")) {
                cli_abort(c("{.pkg birdie} currently supports at most one random intercept term.",
                            "Random slope terms or interacted random effect levels are not supported.",
                            "x"="Found: {.code {re}}."),
                          call=parent.frame())
            }
        }

        if (se_boot > 0) {
            cli_abort(c("Bootstrapping with {.arg se_boot > 0} is only
                        computationally feasible for {.arg model = \"dir\"}.",
                        ">"="If you are using {.arg model = \"auto\"}, you need
                        to remove any random effects terms and specify a completely
                        pooled or fully interacted model to use {.arg model = \"dir\"}."),
                      call=parent.frame())
        }
    } else {
        cli_abort("{.arg model} must be one of {.val {'dir'}},
                  {.val {'mmm'}}, or {.val {'auto'}}.",
                  call=parent.frame())
    }
}

# return TRUE if full interaction structure
check_full_int <- function(tt, covars) {
    int_ord = attr(tt, "order")
    if (length(int_ord) == 0) return(TRUE)
    int_ord = table(factor(int_ord, levels=seq_len(max(int_ord))))
    all(int_ord == choose(length(int_ord), seq_along(int_ord))) # pascal's triangle
}
