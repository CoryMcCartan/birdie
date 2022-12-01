#' @export
birdie <- function(r_probs, formula, data=NULL,
                   prior=NULL, prefix="pr_", ctrl=birdie.ctrl()) {
    # if (missing(data)) cli_abort("{.arg data} must be provided.")

    # figure out type of model and extract response vector
    Y_vec = eval_tidy(f_lhs(formula), data)
    tt = terms(formula)
    covars = all.vars(tt)
    if (!detect_ranef(tt)) {
        d_model = model.frame(formula, data=data, na.action=na.fail)
        if (length(attr(tt, "term.labels")) == 0) { # just Y ~ 1
            method = "pool"
        } else { # saturated
            method = "fixef"
            check_full_int(tt, covars)
        }
    } else {
        d_model = get_all_vars(formula, data=data)
        attr(d_model, "terms") = tt

        method = if (ncol(d_model) == 2) "re1" else "glmm"
    }

    check_covars(r_probs, covars, method) # check predictors
    # set up race probability matrix
    if (!is.matrix(r_probs)) {
        p_rxs = as.matrix(select(r_probs, starts_with(prefix)))
        colnames(p_rxs) = substring(colnames(p_rxs), nchar(prefix)+1L)
    } else {
        p_rxs = r_probs
    }

    # check types
    if (!check_vec(Y_vec)) cli_abort("Response variable must be a character or factor with no missing values.")
    if (nrow(data) != nrow(r_probs))
        cli_abort("{.arg data} and {.arg r_probs} must have the same number of rows.")

    Y_vec = as.factor(Y_vec)
    n_y = nlevels(Y_vec)
    n_r = ncol(p_rxs)

    prior = check_make_prior(prior, method, n_y, n_r)

    # run inference
    t1 <- Sys.time()
    if (method %in% c("pool", "fixef")) {
        res = em_fixef(Y_vec, p_rxs, d_model[-1], prior, ctrl=ctrl)
    } else if (method == "re1") {
        # res = em_re1(Y_vec, p_rxs, d_model[-1], prior, ctrl=ctrl)
        res = em_glmm(Y_vec, p_rxs, d_model, prior, ctrl=ctrl)
        # cli_abort("Method {.val {method}} not yet implemented.")
    } else if (method == "lmer") {
        cli_abort("Method {.val {method}} not yet implemented.")
        # out = em_glmm(Y_vec, r_probs, formula, data, ctrl=ctrl)
    }
    t2 <- Sys.time()

    if (isFALSE(res$converge)) {
        cli_warn(c("EM algorithm did not converge in {ctrl$max_iter} iterations.",
                   ">"="Consider increasing {.arg max_iter}."),
                 call=parent.frame())
    }


    # add names
    colnames(res$map) = colnames(p_rxs)
    rownames(res$map) = levels(Y_vec)

    # format p_ryxs
    colnames(res$p_ryxs) = stringr::str_c(prefix, colnames(p_rxs))
    p_ryxs = as_tibble(res$p_ryxs)
    if (inherits(r_probs, "bisg")) {
        attr(p_ryxs, "S_name") = attr(r_probs, "S_name")
        attr(p_ryxs, "GX_names") = c(names(d_model)[1], attr(r_probs, "GX_names"))
        attr(p_ryxs, "p_r") = attr(r_probs, "p_r")
        attr(p_ryxs, "method") = "birdie"
        class(p_ryxs) = c("bisg", class(p_ryxs))
    }

    # output
    attr(tt, ".Environment") = NULL # save space

    structure(list(
        map = res$map,
        map_sub = res$ests,
        p_ryxs = p_ryxs,
        N = length(Y_vec),
        prior = prior,
        prefix = prefix,
        algo = list(
            method = method,
            iters = res$iters,
            converge = res$converge,
            runtime = as.numeric(t2 - t1, units = "secs")
        ),
        call = match.call()
    ), class="birdie")
}

# Fixed-effects model (includes complete pooling and no pooling)
em_fixef <- function(Y, p_rxs, d_model, prior, ctrl) {
    n_y = nlevels(Y)
    n_r = ncol(p_rxs)
    Y = as.integer(Y)

    # create unique group IDs
    X = to_unique_ids(d_model)
    n_x = max(X)
    est_dim = c(n_r, n_y, n_x)

    # init
    ests = dirichlet_map(Y, X, p_rxs, prior, n_x)
    ests0 <<- ests

    # do EM (accelerated)
    res = ctrl$accel(ests, function(curr) {
        .Call(`_birdie_em_dirichlet`, curr, Y, X, p_rxs, prior, n_x, FALSE)
    }, ctrl, n_x=n_x)

    p_ryxs = calc_bayes(Y, X, res$ests, p_rxs, n_x, n_y)
    ones_mat = matrix(1, nrow=n_y, ncol=n_r)
    est = dirichlet_map(Y, rep_along(Y, 1), p_ryxs, ones_mat, 1) %>%
        matrix(n_y, n_r, byrow=TRUE)

    list(map = est,
         ests = to_array_yrx(res$ests, est_dim),
         iters = res$iters,
         converge = res$converge,
         p_ryxs = p_ryxs)
}

em_re1 <- function(Y, p_rxs, d_model, prior, ctrl) {
    n_y = nlevels(Y)
    n_r = ncol(p_rxs)
    Y = as.integer(Y)
    ones_mat = matrix(1, nrow=n_y, ncol=n_r)

    # create unique group IDs
    X = to_unique_ids(d_model)
    n_x = max(X)
    est_dim = c(n_r, n_y, n_x)

    # init pars and access indices
    p_theta = n_r * n_y * n_x
    p_alpha = n_r * n_y
    idx_theta = seq_len(p_theta)

    idx_alpha = p_theta + seq_len(p_alpha)

    prior_sigma = 0.1
    scale_b = (0.25 - prior_sigma^2) / (2 * prior_sigma^2)

    ests = rep_len(0, p_theta + p_alpha)
    ests[idx_alpha] = log(1.01)

    to_prior = function(ests) {
        matrix(exp(ests[idx_alpha]), n_y, n_r, byrow=TRUE)
    }

    ests[idx_theta] = dirichlet_map(Y, X, p_rxs, to_prior(ests), n_x)

    standata = list(n_y = n_y,
                    n_r = n_r,
                    n_x = n_x,
                    prior_loc_alpha = 2.0,
                    prior_shp_alpha = 100.0)

    # do EM (accelerated)
    res = ctrl$accel(ests, function(curr) {
        # reproject if acceleration has brought us out of bounds
        curr[idx_theta][curr[idx_theta] < 0] = 0 + 1e3*.Machine$double.eps
        curr[idx_theta][curr[idx_theta] > 1] = 1 - 1e3*.Machine$double.eps
        cat(".")

        agg = .Call(`_birdie_em_dirichlet`, curr[idx_theta], Y, X,
                    p_rxs, ones_mat, n_x, TRUE)

        standata$X = to_array_xry(agg, est_dim)

        fit = rstan::optimizing(stanmodels$dir_hier, data=standata,
                                init=list(alpha=t(to_prior(ests))),
                                as_vector=FALSE,
                                tol_obj=0.01, tol_param=ctrl$abstol,
                                tol_rel_obj=10/ctrl$reltol)

        curr[idx_alpha] = log(fit$par$alpha)

        curr[idx_theta] = dirichlet_norm(agg, to_prior(curr), n_x)

        curr
    }, ctrl, n_x=n_x)

    p_ryxs = calc_bayes(Y, X, res$ests[idx_theta], p_rxs, n_x, n_y)
    est = dirichlet_map(Y, rep_along(Y, 1), p_ryxs, ones_mat, 1) %>%
        matrix(n_y, n_r, byrow=TRUE)

    cat('\n')
    print(round(to_prior(res$ests), 2))

    list(map = est,
         ests = to_array_yrx(res$ests[idx_theta], est_dim),
         iters = res$iters,
         converge = res$converge,
         p_ryxs = p_ryxs)

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

em_glmm <- function(Y, p_rxs, d_model, prior, ctrl) {
    # rlang::check_installed("lme4", "to fit a GLMM-based BIRDiE model.")

    n_y = nlevels(Y)
    n_r = ncol(p_rxs)
    Y = as.integer(Y)
    ones_mat = matrix(1, nrow=n_y, ncol=n_r)
    ones = rep_along(Y, 1)

    # create unique group IDs
    X = to_unique_ids(d_model[-1])
    n_x = max(X)
    est_dim = c(n_r, n_y, n_x)

    idx_sub = vctrs::vec_unique_loc(X)
    d_sub = d_model[idx_sub, ][-1]

    ests = dirichlet_map(Y, X, p_rxs, ones_mat * 1.01, n_x) |>
        to_array_xyr(est_dim)

    # form_fit = update.formula(form, cbind(succ, fail) ~ .)
    # form_env = rlang::f_env(form_fit)

    standata = list(
        n_y = n_y,
        N = n_x,
        p = 1L,
        n_grp = n_x,

        X = matrix(1, nrow=n_x),

        grp = X[idx_sub],

        prior_sigma = 0.2,
        prior_beta = 1.0
    )

    res = ctrl$accel(ests, function(curr) {
        cat(".")
        cts = .Call(`_birdie_em_dirichlet`, curr, Y, X,
                    p_rxs, ones_mat, n_x, TRUE) |>
            to_array_xyr(est_dim)
        curr = to_array_xyr(curr, est_dim)

        for (r in seq_len(n_r)) {
            standata$Y = cts[, , r]

            fit = rstan::optimizing(stanmodels$multinom, data=standata,
                                    init=0, check_data=FALSE, as_vector=FALSE,
                                    tol_obj=100*ctrl$abstol, tol_param=ctrl$abstol)
            if (r == 2) fit0 <<- fit

            curr[, , r] = exp(fit$par$lsft)
        }

        to_vec_xyr(curr)
    }, ctrl, n_x=n_x)

    # final global mean
    p_ryxs = calc_bayes(Y, X, res$ests, p_rxs, n_x, n_y)
    est = dirichlet_map(Y, ones, p_ryxs, ones_mat, 1) %>%
        matrix(n_y, n_r, byrow=TRUE)

    list(map = est,
         ests = ests,
         iters = res$iters,
         converge = res$converge,
         p_ryxs = p_ryxs)
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
