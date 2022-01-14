
#' Infer Race Given Name, Location, Outcome, and Covariates
#'
#' Collapsed Gibbs sampler for individual race given name, location, outcome,
#' and covariates.
#'
#' @param X the column containing the outcome
#' @param S the column containing surnames
#' @param G the column containing locations
#' @param W the column(s), if any, containing other covariates. Use `c()` to provide multiple columns
#' @param data the data
#' @param p_rs a data frame, containing a column matching `S` and columns for
#'   each value of `R` giving the conditional probabilities of R given S.
#'   Defaults to a table from the U.S. Census Bureau from 2010.
#' @param p_rgz a data frame, containing columns matching `G`and `W`, and columns for
#'   each value of `R` giving the conditional probabilities of R given W and G.
#'   Defaults to a table from the 2020 decennial census.
#' @param p_r a vector containing the marginal probabilities for each value of
#'   `R`. Defaults to the demographics of North Carolina.
#' @param regularize if `TRUE`, regularize the `S|R` and `G,Z|R` tables
#' @param alpha the number of pseudo-observations in each `X` category
#' @param gibbs if `TRUE`, run the Gibbs sampler
#' @param stan if `TRUE`, run the matrix-based Stan model
#' @param stan_method which inference algorithm to use
#' @param iter the number of Gibbs iterations
#' @param thin how often to record the sampled races
#' @param verbose whether to print additional information
#'
#' @return TBD
#' @export
model_race = function(X, S, G, W=NULL, data=NULL, p_rs=NULL, p_rgz=NULL,
                      p_r=c(white=0.626, black=0.222, hisp=0.098, asian=0.032, other=0.022),
                      regularize=TRUE, alpha=7,
                      gibbs=TRUE, stan=TRUE, stan_method=c("vb", "opt", "hmc"),
                      iter=100, thin=1, verbose=FALSE) {
    check_arg = function(x) {
        x_quo = enquo(x)
        if (!is.numeric(x) || length(x) != 1 || x[1] <= 0)
            cli_abort("{.arg {as_name(x_quo)}} must be a positive number")
    }
    check_arg(alpha)
    check_arg(iter)
    check_arg(thin)

    ## Parse and check input data frame ----------------
    if (missing(data)) cli_abort("{.arg data} must be provided.")
    X_vec = eval_tidy(enquo(X), data)
    S_vec = eval_tidy(enquo(S), data)
    G_vec = eval_tidy(enquo(G), data)
    W_df = data[, eval_select(enquo(W), data)]

    check_vec = function(x) (is.character(x) | is.factor(x)) && !any(is.na(x))

    if (!check_vec(X_vec)) cli_abort("{.arg X} must be a character or factor with no missing values.")
    if (!check_vec(S_vec)) cli_abort("{.arg S} must be a character or factor with no missing values.")
    if (!is.character(G_vec) && !is.factor(G_vec))
        cli_abort("{.arg G} must be a character or factor vector.")
    if (any(is.na(G_vec))) {
        #cli_inform("Missing values found in {.arg G}")
        G_vec = coalesce(G_vec, "<none>")
    }
    if (!all(vapply(W_df, class, character(1)) %in% c("character", "vector"))) {
        cli_abort("{.arg W} must contain only character or factor columns.")
    }
    if (any(is.na(W_df))) cli_abort("Missing values found in {.arg W}")

    X_vec = as.factor(X_vec)
    S_vec = as.factor(S_vec)
    G_vec = as.factor(G_vec)
    GW = cbind(G_vec, W_df)

    if (ncol(W_df) > 0) cli_warn("{.arg W} not yet supported.")

    ## Parse and check input probabilities ----------------
    if (missing(p_rs)) {
        S_vec = as.character(S_vec)
        p_rs = census_surname_table(S_vec, as_name(enquo(S)), p_r, regularize)
        S_vec[!S_vec %in% p_rs[[1]]] = "<generic>"
        S_vec = factor(S_vec, levels=p_rs[[1]])
    } else {
        if (!is.data.frame(p_rs)) cli_abort("{.arg p_rs} must be a data frame.")
        if (ncol(p_rs) != nlevels(R_vec) + 1L)
            cli_abort("Number of racial categories in {.arg p_rs} and {.arg R} must match.")
        if (!all(S_vec %in% p_rs[, 1]))
            cli_abort("Some names are missing from {.arg p_rs}.")
    }
    p_s = prop.table(table(S_vec))
    p_sr = as.matrix(p_rs[, -1])
    for (i in seq_along(p_r)) {
        p_sr[, i] = p_sr[, i] * p_s
        p_sr[, i] = p_sr[, i] / sum(p_sr[, i])
    }

    if (missing(p_rgz)) {
        p_rgz = census_zip_table(G_vec, as_name(enquo(G)), p_r, regularize)
    } else {
        if (!is.data.frame(p_rs)) cli_abort("{.arg p_rgz} must be a data frame.")
        if (!all(colnames(GW) %in% colnames(p_rgz)))
            cli_abort("All columns in {.arg G} and {.arg W} must be in {.arg p_rgz}.")
        if (ncol(p_rgz) != nlevels(R_vec) + ncol(GW))
            cli_abort("Number of racial categories in {.arg p_rgz} and {.arg R} must match.")
        if (nrow(anti_join(GW, p_rgz, by=colnames(GW))) > 0)
            cli_abort("Some {.arg G}/{.arg W} combinations are missing from {.arg p_rgz}.")
    }
    p_g = prop.table(table(G_vec))
    p_gzr = as.matrix(p_rgz[, -1])
    for (i in seq_along(p_r)) {
        p_gzr[, i] = p_gzr[, i] * p_g
        p_gzr[, i] = p_gzr[, i] / sum(p_gzr[, i])
    }

    # TODO remove: debug only
    pp <<- list(
        sr = `rownames<-`(p_sr, levels(S_vec)),
        rs = p_rs,
        gzr = `rownames<-`(p_gzr, levels(G_vec)),
        rgz = p_rgz,
        s = `names<-`(p_s, levels(S_vec)),
        g = `names<-`(p_g, levels(G_vec))
    )


    out = c(est_baseline(X_vec, S_vec, G_vec, p_sr, p_gzr, p_r),
            est_leastsq(X_vec, S_vec, p_sr, p_r))

    if (gibbs) {
        out$gibbs = gibbs_race(iter, thin, X_vec, S_vec, G_vec,
                               log(p_sr), log(p_gzr), log(p_r),
                               matrix(alpha, length(p_r), nlevels(X_vec)),
                               verbose)
    }

    if (stan) {
        out$stan = est_stan(X_vec, S_vec, G_vec, p_sr, p_gzr, p_r,
                            alpha[1], match.arg(stan_method), iter, verbose)
    }

    out
}

est_stan = function(X, S, G, p_sr, p_gzr, p_r, alpha, method, iter, verbose) {
    rlang::check_installed("rstan", "for Stan models")

    stan_data = list(
        N = length(X),
        n_x = nlevels(X),
        n_r = length(p_r),
        n_gz = nlevels(G),
        n_s = nlevels(S),

        X = as.integer(X),
        S = as.integer(S),
        GZ = as.integer(G),

        lp_sr = log(p_sr),
        lp_gzr = log(p_gzr),
        lp_r = log(p_r),
        p_gz = prop.table(table(G)),
        p_x = prop.table(table(X)),

        n_prior_obs = alpha[1]
    )

    if (method == "opt") {
        rstan::optimizing(stanmodels$model_distr, stan_data,
                          verbose=verbose,
                          #draws=iter, importance_resampling=TRUE,
                          tol_grad=1e-6, tol_param=1e-5, init_alpha=1)
    } else if (method == "vb") {
        rstan::vb(stanmodels$model_distr, stan_data,
                  init=0,
                  pars=c("p_xr", "p_r"), algorithm="meanfield",
                  eta=1, adapt_engaged=FALSE,
                  grad_samples=2, eval_elbo=50,
                  iter=700, importance_resampling=FALSE)
    } else {
        rstan::sampling(stanmodels$model_distr, stan_data, pars=c("p_xr", p_r),
                        chains=1, iter=300+iter, warmup=300)
    }
}


est_baseline = function(X, S, G, p_sr, p_gzr, p_r) { # baseline
    m_baseline = matrix(nrow=length(X), ncol=length(p_r))
    for (i in seq_along(p_r)) {
        m_baseline[, i] = p_sr[S, i] * p_gzr[G, i] * p_r[i]
    }
    m_baseline = m_baseline / rowSums(m_baseline)
    colnames(m_baseline) = names(p_r)

    m_baseline_surn = matrix(nrow=length(X), ncol=length(p_r))
    for (i in seq_along(p_r)) {
        m_baseline_surn[, i] = p_sr[S, i] * p_r[i]
    }
    m_baseline_surn = m_baseline_surn / rowSums(m_baseline_surn)
    colnames(m_baseline_surn) = names(p_r)

    list(baseline = m_baseline,
         base_surn = m_baseline_surn)
}

est_leastsq = function(X, S, p_sr, p_r) {
    rlang::check_installed(c("MASS", "glmnet"), "generalized inverse")
    p_rs = p_sr %*% diag(p_r)
    p_rs_inv = t(MASS::ginv(p_rs))
    p_xs = prop.table(table(X, S))

    p_x = prop.table(table(X))
    P_xr_low = outer(p_x, p_r, \(x, y) pmax(0, x+y-1))
    P_xr_high = outer(p_x, p_r, pmin)

    tidy_pred = function(est) {
        est[est > P_xr_high] = P_xr_high[est > P_xr_high]
        est[est < P_xr_low] = P_xr_low[est < P_xr_low]
        rownames(est) = levels(X)
        colnames(est) = names(p_r)
        est
    }

    P_xr_lsq = tidy_pred(p_xs %*% p_rs_inv)
    P_xr_nnls = tidy_pred(rbind(
        as.numeric(glmnet::glmnet(p_rs, t(p_xs[1, ]), lambda=0, lower.limits=0, intercept=FALSE)$beta),
        as.numeric(glmnet::glmnet(p_rs, t(p_xs[2, ]), lambda=0, lower.limits=0, intercept=FALSE)$beta),
        as.numeric(glmnet::glmnet(p_rs, t(p_xs[3, ]), lambda=0, lower.limits=0, intercept=FALSE)$beta),
        as.numeric(glmnet::glmnet(p_rs, t(p_xs[4, ]), lambda=0, lower.limits=0, intercept=FALSE)$beta)
    ) %*% diag(p_r))

    list(lsq = P_xr_lsq,
         nnls = P_xr_nnls)
}


eval_log_score = function(..., eps=1e-9) {
    log_score = function(m, eps=1e-9) {
        R_vec = as.integer(voters$race)
        mean(log(map_dbl(seq_len(nrow(voters)), \(i) m[i, R_vec[i]] + eps)))
    }

    scores = rlang::list2(...) %>%
        vapply(log_score, numeric(1), eps) %>%
        sort(decreasing=TRUE)
    tibble(method = names(scores),
           log_score = scores)
}
