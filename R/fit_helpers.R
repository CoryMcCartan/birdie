
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

# helpers for cat_mixed()
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

# bootstrap helper
weight_maker <- function(N, R, weights) {
    if (N > 1000 && R > 100) {
        function() tabulate(sample.int(N, sum(weights), replace=TRUE), N) / N
    } else { # more computationally intensive but smoother
        function() as.numeric(rdirichlet(1, weights))
    }
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
    terms(update.formula(formula, str_c("~ . -", re_adj)))
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
                       "i"="Missing: {.code {unused}}."),
                     .frequency = "regularly", .frequency_id="birdie_covars_all",
                     call=parent.frame())
        }
    }
}


#' BIRDiE Complete-Data Model Families
#'
#' BIRDiE supports a number of complete-data outcome models, including
#' categorical regression models. Models specific to BIRDiE are listed here.
#' See the Details section of [birdie()] for more information about each
#' model.
#'
#' @param link The link function. Only one option available for categorical
#'   regression models.
#'
#' @returns A list of class `family` containing the specification.

#' @examples
#' cat_dir()
#' cat_mixed()
#'
#' @name birdie-family
NULL

#' @rdname birdie-family
#' @export
cat_dir <- function(link="identity") {
    structure(list(family="cat_dir", link=match.arg(link)),
              class="family")
}
#' @rdname birdie-family
#' @export
cat_mixed <- function(link="softmax") {
    structure(list(family="cat_mixed", link=match.arg(link)),
              class="family")
}


check_model <- function(family, tt, covars, full_int, algorithm) {
    if (family$family == "cat_dir") {
        if (!full_int) {
            x_vars = attr(tt, "term.labels")[attr(tt, "order") == 1]
            cli_abort(c("Fixed effects (no-pooling) model is specified without full
                        interaction structure.",
                        "i"="Models without full interaction structure must be fit
                             using {.arg family = cat_mixed()}.",
                        ">"="For full interaction structure, use
                              `{covars[1]} ~ {paste0(x_vars, collapse=' * ')}`."),
                      call=parent.frame())
        }
        if (count_ranef(tt) > 0) {
            cli_abort("Models with random effects must be fit
                             using {.arg family = cat_mixed()}.",
                      call=parent.frame())
        }
        "cat_dir" # return model type
    } else if (family$family == "cat_mixed") {
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
        if (length(re) == 0 && length(attr(tt, "factors")) == 0) {
            cli_abort("{.fn cat_mixed} requires at least one covariate or random effect.",
                      call=parent.frame())
        }

        if (algorithm == "em_boot") {
            cli_abort("Bootstrapping is only computationally feasible for
                      {.arg family = cat_dir()}.",
                      call=parent.frame())
        }
        if (algorithm == "gibbs") {
            cli_abort("Gibbs sampling is not suppored for {.arg family = cat_dir()}.",
                      call=parent.frame())
        }
        "cat_mixed" # return model type
    } else if (family$family == "gaussian") {
        if (family$link != "identity") {
            cli_abort("Only the identity link is supported for the {.fn gaussian}
                      family.", call=parent.frame())
        }
        "lm"
    } else {
        cli_abort("{.arg family} must be one of {.fn cat_dir} or
                  {.fn cat_mixed}.",
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
