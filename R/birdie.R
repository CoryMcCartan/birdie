#' Fit BIRDiE Models
#'
#' Fits one of two possible Bayesian Instrumental Regression for Disparity
#' Estimation (BIRDiE) models to BISG probabilities and covariates. The simplest
#' Multinomial-Dirichlet model (`dir`) is appropriate when there are no covariates or when
#' all covariates are discrete and fully interacted with another. The more
#' general Multinomial mixed-effects model (`mmm`) is a supports any number of
#' fixed effects and up to one random intercept.
#'
#' `birdie()` uses an expectation-maximization (EM) routine to find the maximum
#' *a posteriori* (MAP) estimate for the specified model. Asymptotic
#' variance-covariance matrices for the MAP estimate are available for the
#' Multinomial-Dirichlet model via bootstrapping (`se_boot`).
#'
#' The Multinomial-Dirichlet model is specified as follows: \deqn{
#'     Y_i \mid R_i, X_i, \Theta \sim \text{Categorical}(\theta_{R_iX_i}) \\
#'     \theta_{rx} \sim \text{Dirichlet}(\alpha_r),
#' } where \eqn{Y} is the outcome variable, \eqn{R} is race, \eqn{X} are
#' covariates (fixed effects), and \eqn{\theta_{rx}} and \eqn{\alpha_r} are
#' vectors with length matching the number of levels of the outcome variable.
#' There is one vector \eqn{\theta_{rx}} for every combination of race and
#' covariates, hence the need for `formula` to either have no covariates or a
#' fully interacted structure.
#'
#' The Multinomial mixed-effects model is specified as follows: \deqn{
#'     Y_i \mid R_i, X_i, \Theta \sim \text{Categorical}(g^{-1}(\mu_{R_iX_i})) \\
#'     \mu_{rxy} = W\beta_{ry} + Zu_{ry} \\
#'     u_{ry} \mid \sigma^2_{ry} \sim \mathcal{N}(0, \sigma^2_{ry}) \\
#'     \beta_{ry} \sim \mathcal{N}(0, s^2_{r\beta}) \\
#'     \sigma_{ry} \sim \text{Gamma}(2, 2/s_{r\sigma}),
#' } where \eqn{\beta_{ry}} are the fixed effects, \eqn{u_{ry}} is the random
#' intercept, and \eqn{g} is a softmax link function.
#' Estimates for \eqn{\beta_{ry}} and \eqn{\sigma_{ry}} are stored in the
#' `$betas` and `$sigmas` elements of the fitted model object.
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
#' @param model A string specifying the type of model to fit: either `"dir"` for
#'   the Multinomial-Dirichlet model or `"mmm"` for the Multinomial
#'   mixed-effects model. The default, `"auto"`, will select the most
#'   computationally efficient model available: `"dir"` if `formula` has no
#'   covariates or a fully-interacted structure, and `"mmm"` otherwise. More
#'   details on the model specifications can be found in the "Details" section
#'   below.
#' @param prior A list with entries specifying the model prior.
#'
#'   When `model="dir"` the only entry is `alpha`, which should be a matrix of
#'   Dirichlet hyperparameters. The matrix should have one row for every level
#'   of the outcome variable and one column for every racial group. The default
#'   prior is a matrix with all entries set to \eqn{1+\epsilon}. When
#'   `model="mmm"`, the `prior` list should contain two scalar entries:
#'   `scale_beta`, the standard deviation on the Normal prior for the fixed
#'   effects, and `scale_sigma`, the prior mean of the standard deviation of the
#'   random intercepts. These can be a single scalar or a vector with an entry
#'   for each racial group.
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
birdie <- function(r_probs, formula, data=NULL, model=c("auto", "dir", "mmm"),
                   prior=NULL, prefix="pr_", se_boot=0, ctrl=birdie.ctrl()) {
    # figure out type of model and extract response vector
    Y_vec = eval_tidy(f_lhs(formula), data)
    tt = terms(formula, keep.order=TRUE)
    covars = all.vars(tt)
    full_int = check_full_int(tt, covars)
    model = match.arg(model)
    if (model == "auto") {
        model = if (count_ranef(tt) == 0 && full_int) "dir" else "mmm"
    }

    # check formula and predictors against model and r_probs
    check_model(model, tt, covars, full_int, se_boot)
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
    if (!check_vec(Y_vec))
        cli_abort("Response variable must be a character or factor with no missing values.")
    if (!is.null(data) && nrow(data) != nrow(r_probs))
        cli_abort("{.arg data} and {.arg r_probs} must have the same number of rows.")

    Y_vec = as.factor(Y_vec)
    n_y = nlevels(Y_vec)
    n_r = ncol(p_rxs)

    prior = check_make_prior(prior, model, levels(Y_vec), races)

    # run inference
    t1 <- Sys.time()
    if (model == "dir") {
        res = em_dir(Y_vec, p_rxs, tt, data, prior, boot=se_boot, ctrl=ctrl)
    } else if (model == "mmm") {
        res = em_mmm(Y_vec, p_rxs, tt, data, prior, ctrl=ctrl)
        # add names
        ex_beta = res$betas[[1]]
        res$betas = array(do.call(cbind, res$betas), dim=c(nrow(ex_beta), n_y, n_r))
        dimnames(res$betas) = list(rownames(ex_beta), levels(Y), races)
        res$betas = aperm(res$betas, c(2L, 3L, 1L))

        names(res$sigmas) = races
        res$sigmas = do.call(cbind, res$sigmas)

        res$linpred = array(do.call(cbind, res$linpred), dim=c(nrow(res$tbl_gx), n_y, n_r))
        dimnames(res$linpred) = list(tbl_gx_names(res$tbl_gx), levels(Y), races)
        res$linpred = aperm(res$linpred, c(2L, 3L, 1L))
    }
    t2 <- Sys.time()

    if (isFALSE(res$converge)) {
        cli_warn(c("EM algorithm did not converge in {ctrl$max_iter} iterations.",
                   ">"="Consider increasing {.arg max_iter} in {.fn birdie.ctrl}."),
                 call=parent.frame())
    }

    # add names
    colnames(res$map) = races
    rownames(res$map) = levels(Y_vec)
    dimnames(res$ests) = c(dimnames(res$map), list(tbl_gx_names(res$tbl_gx)))

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
        N = length(Y_vec),
        betas = if (model == "mmm") res$betas else NULL,
        sigmas = if (model == "mmm") res$sigmas else NULL,
        linpred = if (model == "mmm") res$linpred else NULL,
        prior = prior,
        tbl_gx = as_tibble(res$tbl_gx),
        vec_gx = res$vec_gx,
        y = Y_vec,
        y_name = covars[1],
        prefix = prefix,
        entropy = list(pre = entropy(p_rxs),
                       post = entropy(p_ryxs)),
        algo = list(
            model = model,
            iters = res$iters,
            converge = res$converge,
            runtime = as.numeric(t2 - t1, units = "secs")
        ),
        call = match.call()
    ), class="birdie")
}

# Fixed-effects model (includes complete pooling and no pooling)
em_dir <- function(Y, p_rxs, formula, data, prior, boot, ctrl) {
    d_model = model.frame(formula, data=data, na.action=na.fail)[-1]

    n_y = nlevels(Y)
    n_r = ncol(p_rxs)
    Y = as.integer(Y)

    # create unique group IDs
    X = to_unique_ids(d_model)
    idx_sub = vctrs::vec_unique_loc(X)
    n_x = max(X)
    est_dim = c(n_r, n_y, n_x)

    # init
    ests = dirichlet_map(Y, X, p_rxs, prior$alpha, n_x)

    # do EM (accelerated)
    pb_id = cli::cli_progress_bar("EM iterations", total=NA)
    res = ctrl$accel(ests, function(curr) {
        cli::cli_progress_update(id=pb_id)

        # reproject if acceleration has brought us out of bounds
        curr[curr < 0] = 0 + 1e3*.Machine$double.eps
        curr[curr > 1] = 1 - 1e3*.Machine$double.eps

        .Call(`_birdie_em_dirichlet`, curr, Y, X, p_rxs, prior$alpha, n_x, FALSE)
    }, ctrl, n_x=n_x)
    cli::cli_progress_done(id=pb_id)

    p_ryxs = calc_bayes(Y, X, res$ests, p_rxs, n_x, n_y)
    ones_mat = matrix(1, nrow=n_y, ncol=n_r)
    est = dirichlet_map(Y, rep_along(Y, 1), p_ryxs, ones_mat, 1) %>%
        matrix(n_y, n_r, byrow=TRUE)

    out =  list(map = est,
                ests = to_array_yrx(res$ests, est_dim),
                p_ryxs = p_ryxs,
                tbl_gx = d_model[idx_sub, , drop=FALSE],
                vec_gx = X,
                iters = res$iters,
                converge = res$converge)

    if (boot > 0) {
        boot_ests = boot_dir(res$ests, boot, Y, X, p_rxs, prior, n_x, ctrl)
        out$vcov = cov(t(boot_ests))
    }

    out
}


# Multinomial mixed-effects model
em_mmm <- function(Y, p_rxs, formula, data, prior, ctrl) {
    outcomes = levels(Y)
    n_y = length(outcomes)
    n_r = ncol(p_rxs)
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
    X = model.matrix(fixef_form, data=get_all_vars(fixef_form, data=data))
    N = length(idx_sub)
    if (nrow(X) != length(Y)) {
        cli_abort("Missing values found in data..", call=parent.frame())
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
    ests = dirichlet_map(Y, idx_uniq, p_rxs, ones_mat * 1.0001, n_uniq)
    standata = list(
        n_y = n_y,
        N = N,
        p = ncol(X),
        n_grp = n_grp,

        Y = matrix(0, nrow=nrow(X), ncol=n_y),
        X = X,
        grp = Z,

        prior_sigma = prior$scale_sigma[1],
        prior_beta = prior$scale_beta[1]
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

            sm_ir = get_stanmodel(rstantools_model_multinom, standata)
            fit = optim_model(sm_ir, init=par_l[[r]], skeleton=skeleton,
                              tol_rel_obj=10/ctrl$abstol,
                              tol_obj=10*ctrl$abstol, tol_param=ctrl$abstol)
            all_converged = all_converged && fit$converged
            curr[, r] = fit$par
        }
        last_iter_converge <<- all_converged

        # cat(max(abs(curr - last)), "\n")
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
    est = dirichlet_map(Y, ones, p_ryxs, ones_mat, 1) %>%
        matrix(n_y, n_r, byrow=TRUE)

    out = list(map = est,
               ests = to_array_yrx(ests, est_dim),
               p_ryxs = p_ryxs,
               betas = lapply(par_l, function(x) {
                   out = x$beta
                   colnames(out) = outcomes
                   rownames(out) = colnames(X)
                   out
               }),
               sigmas = lapply(par_l, function(x) setNames(x$sigma_grp, outcomes)),
               linpreds = lapply(par_l, function(x) {
                   m = exp(X %*% x$beta)
                   m / rowSums(m)
               }),
               tbl_gx = d_model[idx_sub, , drop=FALSE],
               vec_gx = idx_uniq,
               iters = res$iters,
               converge = res$converge)
}
