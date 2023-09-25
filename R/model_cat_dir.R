# Fixed-effects model (includes complete pooling and no pooling)
em_cat_dir <- function(Y, p_rxs, formula, data, weights, prior, races, boot, ctrl) {
    d_model = model.frame(formula, data=data, na.action=na.fail)[-1]
    use_w = any(weights != 1.0)

    if (!check_discrete(Y))
        cli_abort("Response variable must be a character or factor with no missing values.",
                  call=parent.frame())
    Y = as.factor(Y)
    nms = levels(Y)
    n_y = nlevels(Y)
    n_r = ncol(p_rxs)
    prior = check_make_prior_cat_dir(prior, Y, p_rxs, "em", races)
    Y = as.integer(Y)

    # create unique group IDs
    X = to_unique_ids(d_model)
    idx_sub = vctrs::vec_unique_loc(X)
    n_x = max(X)
    est_dim = c(n_r, n_y, n_x)

    # init
    ests = dirichlet_map(Y, X, p_rxs * weights, prior$alpha, n_x)

    # do EM (accelerated)
    pb_id = cli::cli_progress_bar("EM iterations", total=NA)
    res = ctrl$accel(ests, function(curr) {
        cli::cli_progress_update(id=pb_id)

        # reproject if acceleration has brought us out of bounds
        curr[curr < 0] = 0 + 1e3*.Machine$double.eps
        curr[curr > 1] = 1 - 1e3*.Machine$double.eps

        if (use_w) {
            .Call(`_birdie_em_dirichlet_wt`, curr, Y, X, weights, p_rxs, prior$alpha, n_x)
        } else {
            .Call(`_birdie_em_dirichlet`, curr, Y, X, p_rxs, prior$alpha, n_x, FALSE)
        }
    }, ctrl, n_x=n_x)
    cli::cli_progress_done(id=pb_id)

    p_ryxs = calc_bayes(Y, X, res$ests, p_rxs, n_x, n_y)
    ones_mat = matrix(1, nrow=n_y, ncol=n_r)
    est = dirichlet_map(Y, rep_along(Y, 1), p_ryxs * weights, ones_mat, 1) %>%
        matrix(n_y, n_r, byrow=TRUE)
    rownames(est) = nms

    out = list(map = est,
               ests = to_array_yrx(res$ests, est_dim),
               p_ryxs = p_ryxs,
               tbl_gx = d_model[idx_sub, , drop=FALSE],
               vec_gx = X,
               prior = prior,
               iters = res$iters,
               converge = res$converge)

    if (boot > 0) {
        boot_ests = boot_cat_dir(res$ests, boot, Y, X, weights, p_rxs, prior, n_x, ctrl)
        out$vcov = cov(t(boot_ests))
    }

    out
}

# bootstrap `em_cat_dir()`
boot_cat_dir <- function(mle, R=10, Y, X, weights, p_rxs, prior, n_x, ctrl) {
    N = length(Y)
    W_tot = sum(weights)
    n_r = ncol(prior$alpha)
    n_y = nrow(prior$alpha)

    out = matrix(nrow=length(prior$alpha), ncol=R)

    ctrl$abstol = max(0.0005, ctrl$abstol)
    ctrl$reltol = max(0.005, ctrl$reltol)
    ctrl$max_iter = 50

    ones = rep_along(Y, 1)
    ones_mat = matrix(1, nrow=n_y, ncol=n_r)
    mk_wt = weight_maker(N, R, weights)

    cli::cli_progress_bar("Bootstrapping", total=R)
    for (i in seq_len(R)) {
        wt = mk_wt()

        res = ctrl$accel(mle, function(curr) {
            # reproject if acceleration has brought us out of bounds
            curr[curr < 0] = 0 + 1e3*.Machine$double.eps
            curr[curr > 1] = 1 - 1e3*.Machine$double.eps

            .Call(`_birdie_em_dirichlet_wt`, curr, Y, X, wt, p_rxs, prior$alpha, n_x)
        }, ctrl, n_x=n_x)

        out[, i] = em_dirichlet_wt(res$ests, Y, ones, wt, p_rxs, ones_mat, 1)

        cli::cli_progress_update()
    }
    cli::cli_progress_done()

    out
}


gibbs_cat_dir <- function(Y, p_rxs, formula, data, weights, prior, races, iter, warmup, ctrl) {
    d_model = model.frame(formula, data=data, na.action=na.fail)[-1]
    use_w = any(weights != 1.0)

    if (!check_discrete(Y))
        cli_abort("Response variable must be a character or factor with no missing values.",
                  call=parent.frame())
    Y = as.factor(Y)
    nms = levels(Y)
    N = length(Y)
    n_y = nlevels(Y)
    n_r = ncol(p_rxs)
    prior = check_make_prior_cat_dir(prior, Y, p_rxs, "gibbs", races)
    alpha_prior_vec = c(prior$alpha)
    Y = as.integer(Y)

    # create unique group IDs
    X = to_unique_ids(d_model)
    idx_sub = vctrs::vec_unique_loc(X)
    ones_vec = rep_along(Y, 1)
    ones_mat = matrix(1, nrow=n_y, ncol=n_r)
    n_x = max(X)
    est_dim = c(n_r, n_y, n_x)

    # init
    if (iter <= 0 || warmup < 0) {
        cli_abort("{.arg iter} and {.arg warmup} must be positive integers.",
                  call=parent.frame())
    }
    ests = matrix(nrow = prod(est_dim), ncol = as.integer(iter + warmup + 1))
    ests_glb = matrix(nrow = n_r*n_y, ncol = ncol(ests))
    p_ryxs = matrix(0, nrow = N, ncol = n_r)
    ests[, 1] = dirichlet_map(Y, X, p_rxs * weights, prior$alpha, n_x)
    n_imp_ish = 50
    imp_ctr = 0
    R_imp = matrix(nrow = N, ncol = n_imp_ish + 1)

    pb_id = cli::cli_progress_bar("Gibbs iterations", total=ncol(ests)-1)
    for (i in seq(2, ncol(ests))) {
        cli::cli_progress_update(id=pb_id)
        p_ryxs_tmp = calc_bayes(Y, X, ests[, i-1], p_rxs, n_x, n_y)
        if (i > 1 + warmup) {
            p_ryxs = p_ryxs + p_ryxs_tmp
            if (i %% (iter / 50) == 0) { # store
                imp_ctr = imp_ctr + 1
                R_imp[, imp_ctr] = mat_rcatp(p_ryxs)
            }
        }
        ests[, i] = gibbs_dir_step(Y, X, weights, p_ryxs_tmp, prior$alpha, n_x)
        ests_glb[, i] = dirichlet_map(Y, ones_vec, p_ryxs_tmp * weights, ones_mat, 1)
    }
    cli::cli_progress_done(id=pb_id)

    p_ryxs = p_ryxs / iter
    idx_use = -seq_len(warmup + 1)
    ests = rowMeans(ests[, idx_use])
    ests_glb = ests_glb[, idx_use]
    est = rowMeans(ests_glb) |>
        matrix(nrow=n_y, ncol=n_r, byrow=TRUE)
    rownames(est) = nms

    list(map = est,
         ests = to_array_yrx(ests, est_dim),
         p_ryxs = p_ryxs,
         vcov = cov(t(ests_glb)),
         tbl_gx = d_model[idx_sub, , drop=FALSE],
         vec_gx = X,
         prior = prior,
         R_imp = R_imp[, seq_len(imp_ctr)],
         iters = iter,
         converge = NA)
}



check_make_prior_cat_dir <- function(prior, Y, p_rxs, algorithm, races) {
    n_r = length(races)
    stopifnot(is.factor(Y))
    n_y = nlevels(Y)
    if (is.null(prior)) {
        cli_inform("Using weakly informative empirical Bayes prior for Pr(Y | R)",
                   .frequency="regularly", .frequency_id="birdie_prior_dir",
                   call=parent.frame())
        ones_mat = matrix(1, nrow=n_y, ncol=n_r)
        est0 = dirichlet_map(Y, rep_along(Y, 1), p_rxs, ones_mat, 1) |>
            matrix(n_y, n_r, byrow=TRUE)
        if (algorithm == "em") {
            prior = list(alpha = ones_mat + est0)
        } else {
            prior = list(alpha = est0)
        }
    } else if (length(prior) == 1 && is.na(prior)) {
        prior = list(
            alpha = matrix(1 + 100*.Machine$double.eps, nrow=n_y, ncol=n_r)
        )
    }

    if (!is.null(colnames(prior$alpha))) {
        prior$alpha = prior$alpha[, match(colnames(prior$alpha), races)]
    }
    if (!is.null(rownames(prior$alpha))) {
        prior$alpha = prior$alpha[match(rownames(prior$alpha), levels(Y)), ]
    }

    if (!"alpha" %in% names(prior) ||
        !is.numeric(prior$alpha) || !is.matrix(prior$alpha) ||
        any(is.na(prior$alpha))) {
        cli_abort(c("With {.arg family=cat_dir()}, {.arg prior} must have an entry
                        {.code alpha} which is a numeric matrix.",
                    "i"="See {.fn birdie::birdie} for details."),
                  call=parent.frame())
    }
    if (nrow(prior$alpha) != n_y) {
        cli_abort("{.arg prior$alpha} must have the same number of rows
                      as there are levels of Y", call=parent.frame())
    }
    if (ncol(prior$alpha) != n_r) {
        cli_abort("{.arg prior$alpha} must have the same number of columns
                      as there are racial groups", call=parent.frame())
    }
    if (any(prior$alpha < 0)) {
        cli_abort("{.arg prior$alpha} must have nonnegative entries", call=parent.frame())
    }
    if (algorithm != "gibbs" && any(prior$alpha <= 1)) {
        cli_warn("A {.arg prior$alpha} with entries that are not
                     strictly greater than 1 may lead to numerical
                     issues.", call=parent.frame())
    }

    prior
}
