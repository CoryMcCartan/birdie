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
        method = if (length(re_terms) == 1) "re1" else "glmm"
    }

    check_covars(r_probs, covars, method) # check predictors
    # set up race probability matrix
    if (!is.matrix(r_probs)) {
        r_probs = as.matrix(select(r_probs, starts_with(prefix)))
        colnames(r_probs) = substring(colnames(r_probs), nchar(prefix)+1L)
    }

    # check types
    if (!check_vec(Y_vec)) cli_abort("{.arg X} must be a character or factor with no missing values.")
    if (nrow(data) != nrow(r_probs))
        cli_abort("{.arg data} and {.arg r_probs} must have the same number of rows.")

    Y_vec = as.factor(Y_vec)
    n_y = nlevels(Y_vec)

    prior = check_make_prior(prior, method, n_y)

    # run inference
    t1 <- Sys.time()
    if (method %in% c("pool", "fixef")) {
        res = em_fixef(Y_vec, r_probs, d_model[-1], prior, ctrl=ctrl)
    } else if (method == "re1") {
        # out = em_re1(Y_vec, r_probs, formula, data, iter=max_iter)
        cli_abort("Method {.val {method}} not yet implemented.")
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
    colnames(res$map) = colnames(r_probs)
    rownames(res$map) = levels(Y_vec)
    colnames(res$p_ryxs) = stringr::str_c(prefix, colnames(r_probs))

    # output
    attr(tt, ".Environment") = NULL # save space

    structure(list(
        map = res$map,
        map_sub = res$ests,
        p_ryxs = as_tibble(res$p_ryxs),
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

    # do EM (accelerated)
    res = SQUAREM::squarem(
        ests,
        function(curr) em_dirichlet(curr, Y, X, p_rxs, prior, n_x),
        control=list(K=3, method="rre", minimize=TRUE, square=TRUE,
                     step.min0=1, step.max0=1, mstep=4, objfn.inc=1,
                     kr=1, tol=ctrl$abstol*(n_x^0.666), maxiter=ctrl$max_iter,
                     trace=FALSE, intermed=FALSE)
    )
    ests = res$par

    p_ryxs = calc_bayes(Y, X, ests, p_rxs, n_x, n_y)
    est = dirichlet_map(Y, rep_along(Y, 1), p_ryxs, rep_len(1, n_y), 1) %>%
        matrix(n_y, n_r, byrow=TRUE)

    list(map = est,
         ests = aperm(array(ests, est_dim), 3:1),
         iters = res$fpevals,
         converge = res$convergence,
         p_ryxs = p_ryxs)
}

em_re1 <- function(Y, p_rxs, form, d_model, prior=rep(1, ncol(p_rxs)),
                   iter=10, tol=0.001) {

}

em_glmm <- function(Y, p_rxs, form, data, ctrl, ref_grp=1L) {
    rlang::check_installed("lme4", "to fit a GLMM-based BIRDiE model.")

    n_y = nlevels(Y)
    n_r = ncol(p_rxs)
    Y = as.integer(Y)
    ones = rep_along(Y, 1)
    # init with (penalized) weighted estimator
    p_ryxs = p_rxs
    est = dirichlet_map(Y, ones, p_rxs, rep(2, n_y), 1)
    p_ryxs = calc_bayes(Y, ones, est, p_rxs, 1, n_y)

    idx_grp = to_unique_ids(lme4::lFormula(form, data=data)$fr[-1])
    n_grp = max(idx_grp)
    grp_revlk = lapply(seq_len(n_grp), \(i) which(idx_grp == i))
    idx_extr = vapply(grp_revlk, function(x) x[1], 1L)
    idx_extr_y = lapply(seq_len(n_y), function(y) {
        intersect(idx_extr, which(Y %in% c(y, ref_grp)))
    })

    ests = array(dim=c(n_grp, n_y, n_r))
    updating = matrix(TRUE, n_y, n_r)

    form_fit = update.formula(form, cbind(succ, fail) ~ .)
    form_env = rlang::f_env(form_fit)

    break_next = FALSE
    cli::cli_progress_bar("Fitting BIRDiE with EM",
                          total=ctrl$max_iter*(n_y-1)*n_r)
    for (i in seq_len(ctrl$max_iter)) {
        last_ests = ests
        # M step
        for (r in seq_len(n_r)) {
            cts = sum_grp(Y, idx_grp, p_ryxs[, r], rep_len(0, n_y), n_y, n_grp)

            for (y in seq_len(n_y)[-ref_grp]) {
                cli::cli_progress_update()
                if (!updating[y, r]) next

                rlang::env_bind(
                    form_env,
                    succ = cts[, y],
                    fail = cts[, ref_grp],
                )
                m = lme4::glmer(form_fit,
                                data = data[idx_extr, ],
                                family=binomial(), nAGQ=0L) %>%
                    suppressWarnings() %>%
                    suppressMessages()

                ests[, y, r] = exp(predict(m, type="link"))
            }

            # normalize
            ests[, ref_grp, r] = 1
            ests[, , r] = ests[, , r] / rowSums(ests[, , r])

            # check convergence
            if (i > 1) {
                for (y in seq_len(n_y)) {
                    if (check_convergence(ests[, y, r], last_ests[, y, r],
                                          ctrl$abstol, ctrl$reltol)) {
                        updating[y, r] = FALSE
                    }
                }
            }
        }

        # E step
        p_ryxs = calc_bayes(Y, idx_grp, ests, p_rxs, n_grp, n_y)

        if (break_next) break
        if (all(!updating)) { # converged -- do one more full update
            updating = matrix(TRUE, n_y, n_r)
            break_next = TRUE
        }
    }
    cli::cli_progress_done()


    # final global mean
    est = matrix(dirichlet_map(Y, ones, p_ryxs, rep_len(1, n_y), 1), n_y, n_r)

    list(map = est,
         ests = ests,
         iters = i,
         converge = break_next,
         p_ryxs = p_ryxs)
}


# Helpers for use in fixed-point iteration accelerators
est_to_vec <- function(ests) {
    as.numeric(do.call(cbind, lapply(ests, as.numeric)))
}
vec_to_ests <- function(vec, n_y, n_r) {
    apply(matrix(vec, nrow=n_y*n_r), 2, function(x) {
        matrix(x, nrow=n_y, ncol=n_r)
    }, simplify=FALSE)
}

#' Control of BIRDiE Model Fitting
#'
#' Constructs control parameters for BIRDiE model fitting.
#' All arguments have defaults.
#'
#' @param max_iter The maximum number of EM iterations.
#' @param abstol The absolute tolerance used in checking convergence.
#' @param reltol The relative tolerance used in checking convergence.
#'
#' @return A list containing the control parameters
#'
#' @examples
#' birdie.ctrl(max_iter=1000)
#'
#' @export
birdie.ctrl <- function(max_iter=1000, abstol=1e-6, reltol=1e-6) {
    stopifnot(max_iter >= 1)
    stopifnot(abstol >= 0)
    stopifnot(reltol >= 0)

    list(max_iter = as.integer(max_iter),
         abstol = abstol,
         reltol = reltol)
}

# Check (and possibly create default) prior
check_make_prior <- function(prior, method, n_y) {
    if (is.null(prior)) {
        prior = rep_len(1, n_y)
        if (method == "pool") {
            cli_inform("Using uniform prior for {.arg alpha} = Pr(X | R)",
                       call=parent.frame())
        } else if (method == "fixef") {
            cli_inform("Using c(1+\u03B5, 1+\u03B5, ... 1+\u03B5) prior for
                       {.arg alpha} = Pr(X | R)", call=parent.frame())
            prior = rep_len(1 + 100*.Machine$double.eps, n_y)
        }
    }
    if (length(prior) != n_y) {
        cli_abort("{.arg alpha} prior must have the same number of elements
                  as there are levels of X", call=parent.frame())
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
        cli_warn(c("Fixed effects (no-pooling) model being fit without
                         full interaction structure. No estimation guarantees.",
                   "i"="For full interaction structure, use
                         `{covars[1]} ~ {paste0(x_vars, collapse=' * ')}`."),
                 call=parent.frame())
    }
}
