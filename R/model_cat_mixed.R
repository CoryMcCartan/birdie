# Categorical mixed-effects model
em_cat_mixed <- function(Y, p_rxs, formula, data, weights, prior, races, ctrl) {
    if (!check_discrete(Y))
        cli_abort("Response variable must be a character or factor with no missing values.",
                  call=parent.frame())
    Y = as.factor(Y)
    nms = levels(Y)
    outcomes = levels(Y)
    n_y = length(outcomes)
    n_r = ncol(p_rxs)
    prior = check_make_prior_cat_mixed(prior, Y, races)
    use_w = any(weights != 1.0)

    Y = as.integer(Y)
    ones_mat = matrix(1, nrow=n_y, ncol=n_r)
    ones = rep_along(Y, 1)

    # find unique rows
    d_model = get_all_vars(formula, data=data)[-1]

    idx_uniq = to_unique_ids(d_model)
    idx_sub = vctrs::vec_unique_loc(idx_uniq)
    n_uniq = max(idx_uniq)
    est_dim = c(n_r, n_y, n_uniq)

    # create fixed effects matrix
    fixef_form = remove_ranef(formula)
    attr(fixef_form, "intercept") = 0 # remove intercept
    X = model.matrix(fixef_form, data=get_all_vars(fixef_form, data=data))
    N = length(idx_sub)
    if (nrow(X) != length(Y)) {
        cli_abort("Missing values found in data.", call=parent.frame())
    }
    X = X[idx_sub, , drop=FALSE]

    # create random effects vector
    if (count_ranef(formula) >= 1) {
        re_expr = attr(formula, "variables")[[2 + which(logi_ranef(formula))]][[3]]
        Z = eval_tidy(re_expr, data=data)
        if (any(is.na(Z))) {
            cli_abort("Missing values found in data.", call=parent.frame())
        }
        Z = to_unique_ids(Z[idx_sub])
        n_grp = max(Z)
    } else {
        Z = rep_along(idx_sub, 1L)
        n_grp = 1L
    }


    # init
    ests = dirichlet_map(Y, idx_uniq, p_rxs * weights, ones_mat * 1.0001, n_uniq)
    standata = list(
        n_y = n_y,
        N = N,
        p = ncol(X),
        n_grp = n_grp,

        Y = matrix(0, nrow=nrow(X), ncol=n_y),
        X = X,
        w = tapply(weights, idx_uniq, sum),
        grp = Z,

        has_int = attr(formula, "intercept"),
        prior_sigma = prior$scale_sigma[1],
        prior_beta = prior$scale_beta[1],
        prior_int = prior$scale_int[1]
    )

    sm = get_stanmodel(rstantools_model_multinom, standata)
    skeleton = get_skeleton(sm)

    n_upar = sm$num_pars_unconstrained()
    par0 = rep_len(0, n_upar*n_r)
    first_iter = TRUE

    pb_id = cli::cli_progress_bar("EM iterations", total=NA)
    last_iter_converge = TRUE
    res = ctrl$accel(par0, function(curr) {
        cli::cli_progress_update(id=pb_id)

        curr = matrix(curr, nrow=n_upar, ncol=n_r)
        par_l = apply(curr, 2, function(x) constrain_pars(sm, skeleton, x))
        ests_vec = if (first_iter) {
            first_iter <<- FALSE
            ests
        } else {
            to_ests_vec(par_l, n_y, n_r, N)
        }

        cts = .Call(`_birdie_em_dirichlet`, ests_vec, Y, idx_uniq,
                    p_rxs, ones_mat, n_uniq, TRUE) %>%
            to_array_xyr(est_dim)

        all_converged = TRUE
        for (r in seq_len(n_r)) {
            standata$Y = cts[, , r]
            standata$prior_sigma = prior$scale_sigma[r]
            standata$prior_beta = prior$scale_beta[r]
            standata$prior_int = prior$scale_int[r]

            sm_ir = get_stanmodel(rstantools_model_multinom, standata)
            fit = optim_model(sm_ir, init=par_l[[r]], skeleton=skeleton,
                              tol_rel_obj=10/ctrl$abstol,
                              tol_obj=10*ctrl$abstol, tol_param=ctrl$abstol)
            all_converged = all_converged && fit$converged
            curr[, r] = fit$par
        }
        last_iter_converge <<- all_converged

        as.numeric(curr)
    }, ctrl, n_x=n_upar*n_r*(4^2)) # extra factor since upar scale is different
    cli::cli_progress_done(id=pb_id)

    if (!last_iter_converge) {
        cli::cli_warn("Final M step did not converge.",
                      call=parent.frame())
    }

    par_l = matrix(res$ests, nrow=n_upar, ncol=n_r) %>%
        apply(2, function(x) constrain_pars(sm, skeleton, x))
    ests = to_ests_vec(par_l, n_y, n_r, N)

    # final global mean and R|YXS
    p_ryxs = calc_bayes(Y, idx_uniq, ests, p_rxs, n_uniq, n_y)
    est = dirichlet_map(Y, ones, p_ryxs * weights, ones_mat, 1) %>%
        matrix(n_y, n_r, byrow=TRUE)
    rownames(est) = nms

    out = list(map = est,
               ests = to_array_yrx(ests, est_dim),
               p_ryxs = p_ryxs,
               beta = lapply(par_l, function(x) {
                   out = x$beta
                   colnames(out) = outcomes
                   rownames(out) = colnames(X)
                   if (standata$has_int == 1) {
                       out = rbind(intercept=x$intercept, out)
                   }
                   out
               }),
               sigma = lapply(par_l, function(x) setNames(x$sigma_grp, outcomes)),
               linpreds = lapply(par_l, function(x) {
                   m = exp(standata$has_int * x$intercept + X %*% x$beta)
                   m / rowSums(m)
               }),
               tbl_gx = d_model[idx_sub, , drop=FALSE],
               vec_gx = NULL,
               prior = prior,
               iters = res$iters,
               converge = res$converge)
}

check_make_prior_cat_mixed <- function(prior, Y, races) {
    n_r = length(races)
    if (is.null(prior)) {
        prior = list(
            scale_int = rep(2, n_r),
            scale_beta = rep(0.2, n_r),
            scale_sigma = rep(0.05, n_r)
        )

        cli_inform(c("Using default prior for Pr(Y | R):",
                     ">"="Prior scale on intercepts:
                              {format(prior$scale_int[1], nsmall=1)}",
                     ">"="Prior scale on fixed effects coefficients:
                              {format(prior$scale_beta[1], nsmall=1)}",
                     ">"="Prior mean of random effects standard deviation:
                              {format(prior$scale_sigma[1], nsmall=2)}"),
                   .frequency="regularly", .frequency_id="birdie_prior_mmm",
                   call=parent.frame())
    }

    if (!is.null(names(prior$scale_beta))) {
        prior$scale_beta = prior$scale_beta[match(names(prior$scale_beta), races)]
    }
    if (!is.null(names(prior$scale_int))) {
        prior$scale_int = prior$scale_int[match(names(prior$scale_int), races)]
    }
    if (!is.null(names(prior$scale_sigma))) {
        prior$scale_sigma = prior$scale_sigma[match(names(prior$scale_sigma), races)]
    }


    if (!all(c("scale_sigma", "scale_int", "scale_beta") %in% names(prior)) ||
        !is.numeric(prior$scale_beta) || !length(prior$scale_beta) %in% c(1L, n_r) ||
        !is.numeric(prior$scale_int) || !length(prior$scale_int) %in% c(1L, n_r) ||
        !is.numeric(prior$scale_sigma) || !length(prior$scale_sigma) %in% c(1L, n_r)) {
        cli_abort(c("With {.arg family=cat_mixed()}, {.arg prior} must have three
                        scalar or vector entries {.code scale_int},
                        {.code scale_beta}, and {.code scale_sigma}.",
                    "i"="See {.fn birdie::birdie} for details."),
                  call=parent.frame())
    }

    if (length(prior$scale_beta) != n_r) {
        prior$scale_beta = rep(prior$scale_beta, n_r)
    }
    if (length(prior$scale_int) != n_r) {
        prior$scale_int = rep(prior$scale_int, n_r)
    }
    if (length(prior$scale_sigma) != n_r) {
        prior$scale_sigma = rep(prior$scale_sigma, n_r)
    }

    prior
}
