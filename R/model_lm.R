# linear model
em_lm <- function(Y, p_rxs, formula, data, weights, prior, races, boot, ctrl) {
    d_model = model.frame(formula, data=data, na.action=na.fail)[-1]
    use_w = any(weights != 1.0)

    if (!(is.numeric(Y) || is.logical(Y)) || any(is.na(Y)))
        cli_abort("Response variable must be numeric with no missing values.",
                  call=parent.frame())
    Y = as.numeric(Y)
    if (vctrs::vec_unique_count(Y) <= max(0.05 * length(Y), 10)) {
        cli_warn(c("Found many duplicate values of the outcome variable.",
                   "i"="A Normal linear model may not be appropriate."),
                 call=parent.frame())
    }
    n_r = ncol(p_rxs)
    prior = check_make_prior_lm(prior, Y, races)

    # model matrix and apply prior
    X = model.matrix(formula, data=get_all_vars(formula, data=data))
    if (nrow(X) != length(Y)) {
        cli_abort("Missing values found in data.", call=parent.frame())
    }
    p = ncol(X)
    Y = c(rep(0, p), Y)
    weights = c(rep(1, p), weights)
    # p_rxs = rbind(matrix(1, nrow=p, ncol=ncol(p_rxs)), p_rxs)
    if (attr(formula, "intercept") == 1) { # has intercept
        if (p > 1) {
            Xsc = scale(X[, -1], center=TRUE, scale=FALSE)
            X[, -1] = Xsc
            X = rbind(diag(c(1/prior$scale_int, rep(1/prior$scale_beta, p - 1))), X)
        } else {
            X = rbind(1/prior$scale_int, X)
        }
    } else {
        Xsc = scale(X, center=TRUE, scale=FALSE)
        X = rbind(diag(rep(1/prior$scale_beta, p)), Xsc)
    }

    # initial estimates
    par0 = lm_mstep(X, Y, p_rxs, weights, use_w, p, prior, ctrl)

    pb_id = cli::cli_progress_bar("EM iterations", total=NA)
    res = ctrl$accel(par0, function(curr) {
        cli::cli_progress_update(id=pb_id)

        coefs = matrix(curr[-1], nrow=p, ncol=n_r)
        p_ryxs = lm_estep(X, Y, coefs, curr[1], p_rxs, p)
        lm_mstep(X, Y, p_ryxs, weights, use_w, p, prior, ctrl)
    }, ctrl, n_x=p)
    cli::cli_progress_done(id=pb_id)

    ests = matrix(res$ests[-1], nrow=p, ncol=n_r)
    est_sigma = res$ests[1]
    p_ryxs = lm_estep(X, Y, ests, est_sigma, p_rxs, p)
    ign = -seq_len(p)
    est = colSums(p_ryxs * (X[ign, ] %*% ests)) / colSums(p_ryxs)

    out = list(map = matrix(est, nrow=1),
               ests = matrix(est, nrow=1),
               p_ryxs = p_ryxs,
               beta = ests,
               sigma = est_sigma,
               linpreds = X[ign, ] %*% ests,
               tbl_gx = X[ign, ],
               vec_gx = seq_len(nrow(X) - p),
               prior = prior,
               iters = res$iters,
               converge = res$converge)

    if (boot > 0) {
        boot_ests = boot_lm(res$ests, boot, Y, X, weights, p_rxs, prior, ctrl)
        out$vcov = cov(t(boot_ests))
    }

    out
}

# helper functions for M and E step for linear model
lm_mstep <- function(X, Y, p_ryxs, weights, use_w, p, prior, ctrl) {
    n_r = ncol(p_ryxs)
    ign = -seq_len(p)
    pars = matrix(nrow=p, ncol=n_r)
    alpha_post = prior$n_sigma + nrow(X) - ncol(X) + 1
    beta_post = prior$loc_sigma^2 * prior$n_sigma

    for (i in seq_len(n_r)) {
        pr = c(rep(1, p), p_ryxs[, i])
        if (use_w) {
            pr = pr * weights
        }
        res = lm.wfit(X, Y, pr, tol=ctrl$abstol)

        pars[, i] = res$coefficient
        beta_post = beta_post + sum(res$residuals[ign]^2 * pr[ign])
    }

    sigma = sqrt(beta_post / alpha_post)
    c(sigma, pars)
}
lm_estep <- function(X, Y, coefs, sigma, p_rxs, p) {
    n_r = ncol(coefs)
    ign = -seq_len(p)
    resid = Y[ign] - X[ign, ] %*% coefs
    p_ryxs = log(p_rxs)
    for (i in seq_len(n_r)) {
        p_ryxs[, i] = p_ryxs[, i] + dnorm(resid[, i], sd=sigma, log=TRUE)
    }
    p_ryxs = exp(safeexpoffset(p_ryxs))
    p_ryxs / rowSums(p_ryxs)
}


boot_lm <- function(mle, R=10, Y, X, weights, p_rxs, prior, ctrl) {
    N = length(Y)
    p = ncol(X)
    ign = -seq_len(p)
    n_r = ncol(p_rxs)
    W_tot = sum(weights)

    out = matrix(nrow=n_r, ncol=R)

    ctrl$abstol = max(0.0005, ctrl$abstol)
    ctrl$reltol = max(0.005, ctrl$reltol)
    ctrl$max_iter = 50

    mk_wt = weight_maker(N, R, weights)

    cli::cli_progress_bar("Bootstrapping", total=R)
    for (i in seq_len(R)) {
        wt = mk_wt()

        res = ctrl$accel(mle, function(curr) {
            coefs = matrix(curr[-1], nrow=p, ncol=n_r)
            p_ryxs = lm_estep(X, Y, coefs, curr[1], p_rxs, p)
            lm_mstep(X, Y, p_ryxs, wt, TRUE, p, prior, ctrl)
        }, ctrl, n_x=p)

        ests = matrix(res$ests[-1], nrow=p, ncol=n_r)
        p_ryxs = lm_estep(X, Y, ests, res$ests[1], p_rxs, p)
        out[, i] = colSums(p_ryxs * (X[ign, ] %*% ests)) / colSums(p_ryxs)

        cli::cli_progress_update()
    }
    cli::cli_progress_done()

    out
}

check_make_prior_lm <- function(prior, Y, races) {
    n_r = length(races)
    if (is.null(prior)) {
        cli_inform("Using weakly informative empirical Bayes prior for Pr(Y | R)",
                   .frequency="regularly", .frequency_id="birdie_prior_dir",
                   call=parent.frame())
        sd_Y = sd(Y)
        prior = list(
            scale_beta = 2.5,
            scale_int = 5 * mean(Y) / sd_Y,
            n_sigma = 5,
            loc_sigma = sd_Y
        )
    }

    if (!all(c("scale_int", "scale_beta", "n_sigma", "loc_sigma") %in% names(prior)) ||
        !is.numeric(prior$scale_beta) || length(prior$scale_beta) != 1 ||
        !is.numeric(prior$scale_int) || length(prior$scale_int) != 1 ||
        !is.numeric(prior$n_sigma) || length(prior$n_sigma) != 1 ||
        !is.numeric(prior$loc_sigma) || length(prior$loc_sigma) != 1) {
        cli_abort(c("With {.arg family=gaussian()}, {.arg prior} must have four
                    scalar or vector entries {.code scale_int}, {.code scale_beta},
                    {.code n_sigma}, and {.code loc_sigma}.",
                    "i"="See {.fn birdie::birdie} for details."),
                  call=parent.frame())
    }

    prior
}
