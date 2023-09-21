#' Fit BIRDiE Models
#'
#' Fits one of two possible Bayesian Instrumental Regression for Disparity
#' Estimation (BIRDiE) models to BISG probabilities and covariates. The simplest
#' Categorical-Dirichlet model ([cat_dir()]) is appropriate when there are no
#' covariates or when all covariates are discrete and fully interacted with
#' another. The more general Categorical mixed-effects model ([cat_mixed()]) is
#' a supports any number of fixed effects and up to one random intercept. For
#' continuous outcomes a Normal linear model is available ([gaussian()]).
#'
#' `birdie()` uses an expectation-maximization (EM) routine to find the maximum
#' *a posteriori* (MAP) estimate for the specified model. Asymptotic
#' variance-covariance matrices for the MAP estimate are available for the
#' Categorical-Dirichlet model via bootstrapping (`se_boot`).
#'
#' The Categorical-Dirichlet model is specified as follows: \deqn{
#'     Y_i \mid R_i, X_i, \Theta \sim \text{Categorical}(\theta_{R_iX_i}) \\
#'     \theta_{rx} \sim \text{Dirichlet}(\alpha_r),
#' } where \eqn{Y} is the outcome variable, \eqn{R} is race, \eqn{X} are
#' covariates (fixed effects), and \eqn{\theta_{rx}} and \eqn{\alpha_r} are
#' vectors with length matching the number of levels of the outcome variable.
#' There is one vector \eqn{\theta_{rx}} for every combination of race and
#' covariates, hence the need for `formula` to either have no covariates or a
#' fully interacted structure.
#'
#' The Categorical mixed-effects model is specified as follows: \deqn{
#'     Y_i \mid R_i, X_i, \Theta \sim \text{Categorical}(g^{-1}(\mu_{R_iX_i})) \\
#'     \mu_{rxy} = W\beta_{ry} + Zu_{ry} \\
#'     u_{ry} \mid \sigma^2_{ry} \sim \mathcal{N}(0, \sigma^2_{ry}) \\
#'     \beta_{ry} \sim \mathcal{N}(0, s^2_{r\beta}) \\
#'     \sigma_{ry} \sim \text{Gamma}(2, 2/s_{r\sigma}),
#' } where \eqn{\beta_{ry}} are the fixed effects, \eqn{u_{ry}} is the random
#' intercept, and \eqn{g} is a softmax link function.
#' Estimates for \eqn{\beta_{ry}} and \eqn{\sigma_{ry}} are stored in the
#' `$beta` and `$sigma` elements of the fitted model object.
#'
#' The Normal linear model is specified as follows: \deqn{
#'     Y_i \mid R_i, \vec X_i, \Theta \sim \mathcal{N}(\vec X_i^\top\vec\theta, \sigma^2) \\
#'     \pi(\sigma^2) \propto \sigma^{-2} \\
#'     \beta \sim \mathcal{N}(0, \Sigma), \\
#' } where \eqn{\vec\theta} is a vector of linear model coefficients.
#' Estimates for \eqn{\theta} and \eqn{\sigma} are stored in the
#' `$beta` and `$sigma` elements of the fitted model object.
#'
#' More details on the models and their properties may be found in the paper
#' referenced below.
#'
#' @param r_probs A data frame or matrix of BISG probabilities, with one row per
#'   individual. The output of [bisg()] can be used directly here.
#' @param formula A two-sided formula object describing the model structure. The
#'   left-hand side is the outcome variable, which must be discrete. A single
#'   random intercept term, denoted with a vertical bar (`"(1 | <term>)"`), is
#'   supported on the right-hand side.
#' @param data An optional data frame containing the variables named in `formula`.
#' @param weights An optional numeric vector specifying likelihood weights.
#' @param family A description of the complete-data model type to fit. Options
#'   are:
#'
#'   - [cat_dir()]: Categorical-Dirichlet model. All covariates must be fully
#'   interacted.
#'   - [cat_mixed()]: Categorical mixed-effects model. Up to one random effect
#'   is supported.
#'   - [gaussian()]: Linear model.
#'
#' See the Details section below for more information on the various models.
#' @param model **Deprecated**.
#'
#' @param prior A list with entries specifying the model prior.
#'
#'   - For the `cat_dir` model, the only entry is `alpha`, which should be a matrix
#'   of Dirichlet hyperparameters. The matrix should have one row for every
#'   level of the outcome variable and one column for every racial group. The
#'   default prior (used when `prior=NULL`) is an empirical Bayes prior equal to
#'   1 plus the weighted-mean estimate of the outcome-race table. A fully
#'   noninformative prior with all entries set to \eqn{1+\epsilon} can be
#'   obtained by setting `prior=NA`.
#'   - For the `cat_mixed` model, the `prior` list should contain three scalar entries:
#'   `scale_int`, the standard deviation on the Normal prior for the intercepts
#'   (which control the global estimates of `Y|R`), `scale_beta`, the standard
#'   deviation on the Normal prior for the fixed effects, and `scale_sigma`, the
#'   prior mean of the standard deviation of the random intercepts. These can be
#'   a single scalar or a vector with an entry for each racial group.
#'   - For the `gaussian` model, the `prior` list  should contain two entries:
#'   `scale_int`, controlling, the standard deviation on the Normal prior for
#'   the intercepts (which control the global estimates of `Y|R`), and
#'   `scale_beta`, controlling the standard deviation on the Normal prior for
#'   the fixed effects. These must be a single scalar. Each is expressed in
#'   terms of the estimated residual standard deviation (i.e., they are
#'   multiplied together to form the "true" prior).
#'
#'   The prior is stored after model fitting in the `$prior` element of the
#'   fitted model object.
#' @param prefix If `r_probs` is a data frame, the columns containing racial
#'   probabilities will be selected as those with names starting with `prefix`.
#'   The default will work with the output of [bisg()].
#' @param se_boot The number of bootstrap replicates to use to compute
#'   approximate standard errors for the main model estimates. Only available
#'   when `model="dir"`. When there are fewer than 1,000 individuals or 100 or
#'   fewer replicates, a Bayesian bootstrap is used instead (i.e., weights are
#'   drawn from a \eqn{\text{Dirichlet}(1, 1, ..., 1)} distribution, which produces more
#'   reliable estimates.
#' @param ctrl A list containing control parameters for the EM algorithm and
#'   optimization routines. A list in the proper format can be made using
#'   [birdie.ctrl()].
#'
#' @return An object of class [`birdie`][birdie::birdie-class], for which many
#'   methods are available. The model estimates may be accessed with
#'   [coef.birdie()], and updated BISG probabilities (conditioning on the
#'   outcome) may be accessed with [fitted.birdie()].
#'
#' @references
#' McCartan, C., Goldin, J., Ho, D.E., & Imai, K. (2022).
#' Estimating Racial Disparities when Race is Not Observed.
#' Available at \url{https://arxiv.org/abs/2303.02580}.
#'
#' @examples
#' data(pseudo_vf)
#'
#' r_probs = bisg(~ nm(last_name) + zip(zip), data=pseudo_vf)
#'
#' # Process zip codes to remove missing values
#' pseudo_vf$zip = proc_zip(pseudo_vf$zip)
#'
#' birdie(r_probs, turnout ~ 1, data=pseudo_vf)
#'
#' birdie(r_probs, turnout ~ zip, data=pseudo_vf)
#'
#' fit = birdie(r_probs, turnout ~ (1 | zip), data=pseudo_vf,
#'        ctrl=birdie.ctrl(abstol=1e-3))
#'
#' summary(fit)
#' coef(fit)
#' fitted(fit)
#'
#' @concept estimators
#' @export
birdie <- function(r_probs, formula, data=NULL, weights=NULL,
                   family=cat_dir(), model=NULL,
                   prior=NULL, prefix="pr_", se_boot=0, ctrl=birdie.ctrl()) {
    # figure out type of model and extract response vector
    Y_vec = eval_tidy(f_lhs(formula), data)
    tt = terms(formula, keep.order=TRUE)
    covars = all.vars(tt)
    full_int = check_full_int(tt, covars)
    # deprecated model arg
    if (!missing(model)) {
        cli_warn(c("{.arg model} was deprecated in {.pkg birdie} 0.4.0.",
                   ">"="Please specify the desired model with the {.arg family} argument instead."))
        family = if (model == "mmm") cat_mixed() else cat_dir()
    }

    # check formula and predictors against model and r_probs
    model = check_model(family, tt, covars, full_int, se_boot)
    check_covars(r_probs, covars, model)

    # set up race probability matrix
    if (is.matrix(r_probs)) {
        p_rxs = r_probs
    } else if (is.data.frame(r_probs)) {
        p_rxs = as.matrix(select(r_probs, starts_with(prefix)))
    } else {
        cli_abort("{.arg r_probs} must be a matrix or data frame.")
    }
    races = stringr::str_sub(colnames(p_rxs), nchar(prefix)+1L)

    # check types
    if (!is.null(data) && nrow(data) != nrow(r_probs))
        cli_abort("{.arg data} and {.arg r_probs} must have the same number of rows.")

    n_r = ncol(p_rxs)
    if (is.null(weights)) {
        weights = rep_along(Y_vec, 1.0)
    } else if (!is.numeric(weights) && length(weights) != length(Y_vec)) {
        cli_abort("{.arg weights} must be a numeric vector with one entry for each observation.")
    }

    # run inference
    t1 <- Sys.time()
    if (model == "cat_dir") {
        res = em_cat_dir(Y_vec, p_rxs, tt, data, weights,
                         prior, races, boot=se_boot, ctrl=ctrl)
        # add names
        dimnames(res$ests) = c(dimnames(res$map), list(tbl_gx_names(res$tbl_gx)))
    } else if (model == "cat_mixed") {
        res = em_cat_mixed(Y_vec, p_rxs, tt, data, weights,
                           prior, races, ctrl=ctrl)
        # add names
        n_y = nrow(res$map)
        ex_beta = res$beta[[1]]
        res$beta = array(do.call(cbind, res$beta), dim=c(nrow(ex_beta), n_y, n_r))
        dimnames(res$beta) = list(rownames(ex_beta), levels(Y_vec), races)
        res$beta = aperm(res$beta, c(2L, 3L, 1L))

        names(res$sigma) = races
        res$sigma = do.call(cbind, res$sigma)

        res$linpred = array(do.call(cbind, res$linpred), dim=c(nrow(res$tbl_gx), n_y, n_r))
        dimnames(res$linpred) = list(tbl_gx_names(res$tbl_gx), levels(Y_vec), races)
        res$linpred = aperm(res$linpred, c(2L, 3L, 1L))

        dimnames(res$ests) = c(dimnames(res$map), list(tbl_gx_names(res$tbl_gx)))
    } else if (model == "lm") {
        res = em_lm(Y_vec, p_rxs, tt, data, weights,
                    prior, races, boot=se_boot, ctrl=ctrl)
        # add names
        rownames(res$map) = covars[1]
    }
    t2 <- Sys.time()

    if (isFALSE(res$converge)) {
        cli_warn(c("EM algorithm did not converge in {ctrl$max_iter} iterations.",
                   ">"="Consider increasing {.arg max_iter} in {.fn birdie.ctrl}."),
                 call=parent.frame())
    }

    # format outputs, p_ryxs
    colnames(res$map) = races
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
        N = length(Y_vec),
        beta = if (model %in% c("cat_mixed", "lm")) res$beta else NULL,
        sigma = if (model %in% c("cat_mixed", "lm")) res$sigma else NULL,
        linpred = if (model %in% c("cat_mixed", "lm")) res$linpred else NULL,
        prior = res$prior,
        tbl_gx = as_tibble(res$tbl_gx),
        vec_gx = res$vec_gx,
        y = Y_vec,
        y_name = covars[1],
        prefix = prefix,
        entropy = list(pre = entropy(p_rxs),
                       post = entropy(p_ryxs)),
        algo = list(
            model = model,
            family = family$family,
            iters = res$iters,
            converge = res$converge,
            runtime = as.numeric(t2 - t1, units = "secs"),
            version = as.character(packageVersion("birdie"))
        ),
        call = match.call()
    ), class="birdie")
}

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
    prior = check_make_prior_cat_dir(prior, Y, p_rxs, races)
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
