#' @export
birdie <- function(r_probs, formula, data=NULL,
                   prior=NULL, prefix="pr_", se_boot=0, ctrl=birdie.ctrl()) {
    # figure out type of model and extract response vector
    Y_vec = eval_tidy(f_lhs(formula), data)
    tt = terms(formula)
    covars = all.vars(tt)
    if (!detect_ranef(tt)) {
        d_model = model.frame(formula, data=data, na.action=na.fail)
        method = "fixef"
        if (length(attr(tt, "term.labels")) > 0) { # not just Y ~ 1
            check_full_int(tt, covars)
        }
    } else {
        d_model = get_all_vars(formula, data=data)
        attr(d_model, "terms") = tt
        method = "mmm"
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
    if (!check_vec(Y_vec))
        cli_abort("Response variable must be a character or factor with no missing values.")
    if (!is.null(data) && nrow(data) != nrow(r_probs))
        cli_abort("{.arg data} and {.arg r_probs} must have the same number of rows.")

    Y_vec = as.factor(Y_vec)
    n_y = nlevels(Y_vec)
    n_r = ncol(p_rxs)

    prior = check_make_prior(prior, method, n_y, n_r)

    # run inference
    t1 <- Sys.time()
    if (method %in% c("pool", "fixef")) {
        res = em_fixef(Y_vec, p_rxs, d_model[-1], prior, boot=se_boot, ctrl=ctrl)
    } else if (method == "mmm") {
        res = em_mmm(Y_vec, p_rxs, d_model, prior, ctrl=ctrl)
    }
    t2 <- Sys.time()

    if (isFALSE(res$converge)) {
        cli_warn(c("EM algorithm did not converge in {ctrl$max_iter} iterations.",
                   ">"="Consider increasing {.arg max_iter} in {.fn birdie.ctrl}."),
                 call=parent.frame())
    }

    # add names
    colnames(res$map) = colnames(p_rxs)
    rownames(res$map) = levels(Y_vec)

    # format p_ryxs
    colnames(res$p_ryxs) = stringr::str_c(prefix, colnames(p_rxs))
    p_ryxs = as_tibble(res$p_ryxs)
    if (inherits(r_probs, "bisg")) {
        p_ryxs = reconstruct.bisg(p_ryxs, r_probs, cat_nms=names(d_model)[1])
    }

    # output
    attr(tt, ".Environment") = NULL # save space

    structure(list(
        map = res$map,
        map_sub = res$ests,
        p_ryxs = p_ryxs,
        vcov = if (se_boot > 0) res$vcov else NULL,
        se = if (se_boot > 0) vcov_to_se(res$vcov, res$map) else NULL,
        fit = res$fit,
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
em_fixef <- function(Y, p_rxs, d_model, prior, boot, ctrl) {
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
    res = ctrl$accel(ests, function(curr) {
        # reproject if acceleration has brought us out of bounds
        curr[curr < 0] = 0 + 1e3*.Machine$double.eps
        curr[curr > 1] = 1 - 1e3*.Machine$double.eps

        .Call(`_birdie_em_dirichlet`, curr, Y, X, p_rxs, prior, n_x, FALSE)
    }, ctrl, n_x=n_x)

    p_ryxs = calc_bayes(Y, X, res$ests, p_rxs, n_x, n_y)
    ones_mat = matrix(1, nrow=n_y, ncol=n_r)
    est = dirichlet_map(Y, rep_along(Y, 1), p_ryxs, ones_mat, 1) %>%
        matrix(n_y, n_r, byrow=TRUE)

    out =  list(map = est,
                ests = to_array_yrx(res$ests, est_dim),
                iters = res$iters,
                converge = res$converge,
                p_ryxs = p_ryxs)

    if (boot > 0) {
        boot_ests = boot_fixef(res$ests, boot, Y, X, p_rxs, prior, n_x, ctrl)
        out$vcov = cov(t(boot_ests))
    }

    out
}

boot_fixef <- function(mle, R=10, Y, X, p_rxs, prior, n_x, ctrl) {
    N = length(Y)
    n_r = ncol(prior)
    n_y = nrow(prior)

    out = matrix(nrow=length(prior), ncol=R)

    ctrl$abstol = 0.0005 #min(ctrl$abstol * 1000, 0.001)
    ctrl$reltol = 0.005 #min(ctrl$reltol * 1000, 0.001)
    ctrl$max_iter = 50

    ones = rep_along(Y, 1)
    ones_mat = matrix(1, nrow=n_y, ncol=n_r)
    iters = numeric(R)

    if (N > 1000 && R > 100) {
        mk_wt = function() tabulate(sample.int(N, N, replace=TRUE), N) / N
    } else { # more computationally intensive but smoother
        mk_wt = function() as.numeric(rdirichlet(1, ones))
    }

    cli::cli_progress_bar("Bootstrapping", total=R)
    for (i in seq_len(R)) {
        wt = mk_wt()

        res = ctrl$accel(mle, function(curr) {
            # reproject if acceleration has brought us out of bounds
            curr[curr < 0] = 0 + 1e3*.Machine$double.eps
            curr[curr > 1] = 1 - 1e3*.Machine$double.eps

            .Call(`_birdie_em_dirichlet_wt`, curr, Y, X, wt, p_rxs, prior, n_x)
        }, ctrl, n_x=n_x)
        iters[i] = res$iters

        out[, i] = em_dirichlet_wt(res$ests, Y, ones, wt, p_rxs, ones_mat, 1)

        cli::cli_progress_update()
    }
    cli::cli_progress_done()

    out
}


em_mmm <- function(Y, p_rxs, d_model, prior, ctrl) {
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

    out = list(map = est,
               ests = ests,
               iters = res$iters,
               converge = res$converge,
               p_ryxs = p_ryxs)
}
