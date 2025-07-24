#' Fit BIRDiE Models
#'
#' Fits one of three possible Bayesian Instrumental Regression for Disparity
#' Estimation (BIRDiE) models to BISG probabilities and covariates. The simplest
#' Categorical-Dirichlet model ([cat_dir()]) is appropriate when there are no
#' covariates or when all covariates are discrete and fully interacted with
#' another. The more general Categorical mixed-effects model ([cat_mixed()]) is
#' a supports any number of fixed effects and up to one random intercept. For
#' continuous outcomes a Normal linear model is available ([gaussian()]).
#'
#' By default, `birdie()` uses an expectation-maximization (EM) routine to find
#' the maximum *a posteriori* (MAP) estimate for the specified model. Asymptotic
#' variance-covariance matrices for the MAP estimate are available for the
#' Categorical-Dirichlet and Normal linear models via bootstrapping.
#' Full Bayesian inference is supported via Gibbs sampling for the
#' Categorical-Dirichlet and Normal linear models as well.
#'
#' Whatever model or method is used, a finite-population estimate of the
#' outcome-given-race distribution for the entire observed sample is always
#' calculated and stored as `$est` in the returned object, which can be accessed
#' with [coef.birdie()] as well.
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
#'     u_{r} \mid \vec\sigma_{r}, L_r \sim \mathcal{N}(0,
#'      \text{diag}(\vec\sigma_{r})C_r\text{diag}(\vec\sigma_{r})) \\
#'     \beta_{ry} \sim \mathcal{N}(0, s^2_{r\beta}) \\
#'     \sigma_{ry} \sim \text{Inv-Gamma}(4, 3s_{r\sigma}) \\
#'     C_r \sim \text{LKJ}(2),
#' } where \eqn{\beta_{ry}} are the fixed effects, \eqn{u_{ry}} is the random
#' intercept, and \eqn{g} is a softmax link function.
#' Estimates for \eqn{\beta_{ry}} and \eqn{\sigma_{ry}} are stored in the
#' `$beta` and `$sigma` elements of the fitted model object.
#'
#' The Normal linear model is specified as follows: \deqn{
#'     Y_i \mid R_i, \vec X_i, \Theta \sim \mathcal{N}(\vec X_i^\top\vec\theta, \sigma^2) \\
#'     \sigma^2 \sim \text{Inv-Gamma}(n_\sigma/2, l_\sigma^2 n_\sigma/2) \\
#'     \beta_{\text{intercept}} \sim \mathcal{N}(0, s^2_\text{int}) \\
#'     \beta_k \sim \mathcal{N}(0, s^2_\beta), \\
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
#' @param prior A list with entries specifying the model prior.
#'
#'   - For the `cat_dir` model, the only entry is `alpha`, which should be a matrix
#'   of Dirichlet hyperparameters. The matrix should have one row for every
#'   level of the outcome variable and one column for every racial group. The
#'   default prior (used when `prior=NULL`) is an empirical Bayes prior equal to
#'   the weighted-mean estimate of the outcome-race table. A fully
#'   noninformative prior with all entries set to \eqn{\epsilon} can be obtained
#'   by setting `prior=NA`. When `prior=NULL` and `algorithm="em"` or
#'   `"em_boot"`, 1 is added to the prior so that the posterior mode, rather
#'   than the mean, is shrunk toward these values.
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
#' @param algorithm The inference algorithm to use. One of 3 options:
#'
#'   - `"em"`: An expectation-maximization algorithm which will perform inference
#'   for the maximum a posteriori (MAP) parameter values. Computationally
#'   efficient and supported by all the model families. No uncertainty
#'   quantification.
#'   - `"gibbs"`: A Gibbs sampler for performing full Bayesian inference.
#'   Generally more computationally demanding than the EM algorithm, but
#'   provides uncertainty quantification. Currently supported for `cat_dir()` and
#'   `gaussian()` model families. Computation-reliability tradeoff can
#'   be controlled with `iter` argument.
#'   - `"em_boot"`: Bootstrapped version of EM algorithm. Number of bootstrap
#'   replicates controlled by `iter` parameter. Provides approximate uncertainty
#'   quantification. Currently supported for `cat_dir()` and
#'   `gaussian()` model families.
#' @param iter The number of post-warmup Gibbs samples, or the number of
#'   bootstrap replicates to use to compute approximate standard errors for the
#'   main model estimates. Only available when `family=cat_dir()` or
#'   `gaussian()`. Ignored if `algorithm="em"`.
#'
#'   For bootstrapping, when there are fewer than 1,000 individuals or 100 or
#'   fewer replicates, a Bayesian bootstrap is used instead (i.e., weights are
#'   drawn from a \eqn{\text{Dirichlet}(1, 1, ..., 1)} distribution, which
#'   produces more reliable estimates.
#' @param warmup Number of warmup iterations for Gibbs sampling.
#' Ignored unless `algorithm="gibbs"`.
#' @param prefix If `r_probs` is a data frame, the columns containing racial
#'   probabilities will be selected as those with names starting with `prefix`.
#'   The default will work with the output of [bisg()].
#' @param ctrl A list containing control parameters for the EM algorithm and
#'   optimization routines. A list in the proper format can be made using
#'   [birdie.ctrl()].
#'
#' @return An object of class [`birdie`][birdie::birdie-class], for which many
#'   methods are available. The model estimates may be accessed with
#'   [coef.birdie()], and updated BISG probabilities (conditioning on the
#'   outcome) may be accessed with [fitted.birdie()]. Uncertainty estimates, if
#'   available, can be accessed with `$se` and [vcov.birdie()].
#'
#' @references
#' McCartan, C., Fisher, R., Goldin, J., Ho, D.E., & Imai, K. (2025).
#' Estimating Racial Disparities when Race is Not Observed. *Journal of the
#' American Statistical Association*. Available at
#' \doi{10.1080/01621459.2025.2526695}.
#'
#' @examples
#' \donttest{
#' data(pseudo_vf)
#'
#' r_probs = bisg(~ nm(last_name) + zip(zip), data=pseudo_vf)
#'
#' # Process zip codes to remove missing values
#' pseudo_vf$zip = proc_zip(pseudo_vf$zip)
#'
#' fit = birdie(r_probs, turnout ~ 1, data=pseudo_vf)
#' print(fit)
#' fit$se # uncertainty quantification
#'
#' fit = birdie(r_probs, turnout ~ zip, data=pseudo_vf, algorithm="gibbs")
#'
#' fit = birdie(r_probs, turnout ~ (1 | zip), data=pseudo_vf,
#'              family=cat_mixed(), ctrl=birdie.ctrl(abstol=1e-3))
#'
#' summary(fit)
#' coef(fit)
#' fitted(fit)
#' }
#'
#' @concept estimators
#' @export
birdie <- function(r_probs, formula, data, family=cat_dir(), prior=NULL, weights=NULL,
                   algorithm=c("em", "gibbs", "em_boot"), iter=400, warmup=50,
                   prefix="pr_", ctrl=birdie.ctrl()) {
    # figure out type of model and extract response vector
    Y_vec = eval_tidy(f_lhs(formula), data)
    tt = terms(formula, keep.order=TRUE)
    covars = all.vars(tt)
    full_int = check_full_int(tt, covars)

    # check formula and predictors against model and r_probs
    algorithm = match.arg(algorithm)
    model = check_model(family, tt, covars, full_int, algorithm)
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
    if (missing(data)) data = NULL
    if (!is.null(data) && nrow(data) != nrow(r_probs))
        cli_abort("{.arg data} and {.arg r_probs} must have the same number of rows.")

    n_r = ncol(p_rxs)
    weights = check_make_weights(weights, Y_vec)

    # run inference
    t1 <- Sys.time()
    if (model == "cat_dir") {
        if (algorithm == "gibbs") {
            res = gibbs_cat_dir(Y_vec, p_rxs, tt, data, weights,
                                prior, races, iter, warmup, ctrl)
        } else {
            se_boot = if (algorithm == "em_boot") iter else 0L
            res = em_cat_dir(Y_vec, p_rxs, tt, data, weights,
                             prior, races, boot=se_boot, ctrl=ctrl)
        }
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
        if (algorithm == "gibbs") {
            res = gibbs_lm(Y_vec, p_rxs, tt, data, weights,
                        prior, races, iter, warmup, ctrl)
        } else {
            se_boot = if (algorithm == "em_boot") iter else 0L
            res = em_lm(Y_vec, p_rxs, tt, data, weights,
                        prior, races, boot=se_boot, ctrl=ctrl)
        }
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
    dimnames(res$ests)[[2]] = races
    colnames(res$p_ryxs) = colnames(p_rxs)
    p_ryxs = as_tibble(res$p_ryxs)
    if (inherits(r_probs, "bisg")) {
        p_ryxs = reconstruct_bisg(p_ryxs, r_probs, cat_nms=covars[1])
    }

    # output
    attr(tt, ".Environment") = NULL # save space

    structure(list(
        map = res$map,
        map_sub = res$ests,
        p_ryxs = p_ryxs,
        vcov = if (algorithm != "em") res$vcov else NULL,
        se = if (algorithm != "em") vcov_to_se(res$vcov, res$map) else NULL,
        N = length(Y_vec),
        beta = if (model %in% c("cat_mixed", "lm")) res$beta else NULL,
        sigma = if (model %in% c("cat_mixed", "lm")) res$sigma else NULL,
        linpred = if (model %in% c("cat_mixed", "lm")) res$linpred else NULL,
        prior = res$prior,
        tbl_gx = as_tibble(res$tbl_gx),
        vec_gx = res$vec_gx,
        R_imp = if (algorithm == "gibbs") res$R_imp else NULL,
        y = if (model == "lm") Y_vec else as.factor(Y_vec),
        y_name = covars[1],
        prefix = prefix,
        entropy = list(pre = entropy(p_rxs),
                       post = entropy(p_ryxs)),
        algo = list(
            model = model,
            algorithm = algorithm,
            family = family$family,
            iters = if (algorithm != "em") iter else res$iters,
            converge = if (algorithm != "gibbs") res$converge else NULL,
            runtime = as.numeric(t2 - t1, units = "secs"),
            version = as.character(packageVersion("birdie"))
        ),
        call = match.call()
    ), class="birdie")
}
