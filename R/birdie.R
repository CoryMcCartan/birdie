#' @export
birdie <- function(r_probs, formula, data=NULL, method=c("auto", "fixef", "mmm"),
                   prior=NULL, prefix="pr_", se_boot=0, ctrl=birdie.ctrl()) {
    # figure out type of model and extract response vector
    Y_vec = eval_tidy(f_lhs(formula), data)
    tt = terms(formula, keep.order=TRUE)
    covars = all.vars(tt)
    full_int = check_full_int(tt, covars)
    method = match.arg(method)
    if (method == "auto") {
        method = if (count_ranef(tt) == 0 && full_int) "fixef" else "mmm"
    }

    # check formula and predictors against method and r_probs
    check_method(method, tt, covars, full_int, se_boot)
    check_covars(r_probs, covars, method)

    # set up race probability matrix
    if (!is.matrix(r_probs)) {
        p_rxs = as.matrix(select(r_probs, starts_with(prefix)))
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
    if (method == "fixef") {
        res = em_fixef(Y_vec, p_rxs, tt, data, prior, boot=se_boot, ctrl=ctrl)
    } else if (method == "mmm") {
        res = em_mmm(Y_vec, p_rxs, tt, data, prior, ctrl=ctrl)
    }
    t2 <- Sys.time()

    if (isFALSE(res$converge)) {
        cli_warn(c("EM algorithm did not converge in {ctrl$max_iter} iterations.",
                   ">"="Consider increasing {.arg max_iter} in {.fn birdie.ctrl}."),
                 call=parent.frame())
    }

    # add names
    colnames(res$map) = substring(colnames(p_rxs), nchar(prefix)+1L)
    rownames(res$map) = levels(Y_vec)

    # format p_ryxs
    colnames(res$p_ryxs) = colnames(p_rxs)
    p_ryxs = as_tibble(res$p_ryxs)
    if (inherits(r_probs, "bisg")) {
        p_ryxs = reconstruct.bisg(p_ryxs, r_probs, cat_nms=covars[1])
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
em_fixef <- function(Y, p_rxs, formula, data, prior, boot, ctrl) {
    d_model = model.frame(formula, data=data, na.action=na.fail)[-1]

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
    pb_id = cli::cli_progress_bar("EM iterations", total=NA)
    res = ctrl$accel(ests, function(curr) {
        cli::cli_progress_update(id=pb_id)

        # reproject if acceleration has brought us out of bounds
        curr[curr < 0] = 0 + 1e3*.Machine$double.eps
        curr[curr > 1] = 1 - 1e3*.Machine$double.eps

        .Call(`_birdie_em_dirichlet`, curr, Y, X, p_rxs, prior, n_x, FALSE)
    }, ctrl, n_x=n_x)
    cli::cli_progress_done(id=pb_id)

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

# Multinomial mixed-effects model
em_mmm <- function(Y, p_rxs, formula, data, prior, ctrl) {
    n_y = nlevels(Y)
    n_r = ncol(p_rxs)
    Y = as.integer(Y)
    ones_mat = matrix(1, nrow=n_y, ncol=n_r)
    ones = rep_along(Y, 1)

    # find unique rows
    d_model = get_all_vars(formula, data=data)[-1]
    if (any(is.na(d_model))) {
        cli_abort("Missing values found in data.", call=parent.frame())
    }
    idx_uniq = to_unique_ids(d_model)
    idx_sub = vctrs::vec_unique_loc(idx_uniq)
    n_uniq = max(idx_uniq)
    est_dim = c(n_r, n_y, n_uniq)

    # create fixed effects matrix
    fixef_form = remove_ranef(formula)
    X = model.matrix(fixef_form, data=data)[idx_sub, , drop=FALSE]

    # create random effects vector
    if (count_ranef(formula) >= 1) {
        Z = to_unique_ids(d_model[idx_sub, which(logi_ranef(formula))])
        n_grp = max(Z)
    } else {
        Z = rep_along(idx_sub, 1L)
        n_grp = 1L
    }

    # init
    ests = dirichlet_map(Y, idx_uniq, p_rxs, ones_mat * 1.0001, n_uniq) |>
        to_array_ryx(est_dim)
    standata = list(
        n_y = n_y,
        N = nrow(X),
        p = ncol(X),
        n_grp = n_grp,

        Y = matrix(0, nrow=nrow(X), ncol=n_y),
        X = X,
        grp = Z,

        prior_sigma = prior$sigma,
        prior_beta = prior$beta
    )

    sm <<- rstan::sampling(stanmodels$multinom, data=standata, chains=0) |>
        suppressMessages()

    n_upar = rstan::get_num_upars(sm)
    par0 = rep_len(0, n_upar*n_r)

    pb_id = cli::cli_progress_bar("EM iterations", total=NA)
    res = ctrl$accel(par0, function(curr) {
        cli::cli_progress_update(id=pb_id)

        curr = matrix(curr, nrow=n_upar, ncol=n_r)
        par_l = apply(curr, 2, fn_constr(sm))

        cts = .Call(`_birdie_em_dirichlet`, ests, Y, idx_uniq,
                    p_rxs, ones_mat, n_uniq, TRUE) |>
            to_array_xyr(est_dim)
        ests_arr = to_array_xyr(ests, est_dim)

        for (r in seq_len(n_r)) {
            standata$Y = cts[, , r]

            fit = rstan::optimizing(stanmodels$multinom, data=standata,
                                    init=par_l[[r]], check_data=FALSE, as_vector=FALSE,
                                    tol_obj=10*ctrl$abstol, tol_param=ctrl$abstol)
            curr[, r] = rstan::unconstrain_pars(sm, fit$par)
            if (any(is.nan(curr))) browser()

            ests_arr[, , r] = exp(fit$par$lsft)
            if (any(is.nan(ests_arr))) browser()
        }

        ests <<- to_vec_xyr(ests_arr)
        as.numeric(curr)
    }, ctrl, n_x=n_upar)
    cli::cli_progress_done(id=pb_id)

    # final global mean and R|YXS
    p_ryxs = calc_bayes(Y, idx_uniq, ests, p_rxs, n_uniq, n_y)
    est = dirichlet_map(Y, ones, p_ryxs, ones_mat, 1) %>%
        matrix(n_y, n_r, byrow=TRUE)

    out = list(map = est,
               ests = to_array_xyr(ests, est_dim),
               iters = res$iters,
               converge = res$converge,
               p_ryxs = p_ryxs)
}
