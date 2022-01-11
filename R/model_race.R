
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
#' @param p_rwg a data frame, containing columns matching `G`and `W`, and columns for
#'   each value of `R` giving the conditional probabilities of R given W and G.
#'   Defaults to a table from the 2020 decennial census.
#' @param p_r a vector containing the marginal probabilities for each value of
#'   `R`. Defaults to the demographics of North Carolina.
#' @param alpha the number of pseudo-observations in each `X` category
#' @param iter the number of Gibbs iterations
#' @param thin how often to record the sampled races
#' @param verbose whether to print additional information
#'
#' @return TBD
#' @export
model_race = function(X, S, G, W=NULL, data=NULL, p_rs=NULL, p_rwg=NULL,
                      p_r=c(white=0.626, black=0.222, hisp=0.098, asian=0.032, other=0.022),
                      alpha=1,
                      gibbs=TRUE, stan=TRUE,
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
        p_rs = census_surname_table(S_vec, as_name(enquo(S)), p_r)
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

    if (missing(p_rwg)) {
        p_rwg = census_zip_table(G_vec, as_name(enquo(G)), p_r)
    } else {
        if (!is.data.frame(p_rs)) cli_abort("{.arg p_rwg} must be a data frame.")
        if (!all(colnames(GW) %in% colnames(p_rwg)))
            cli_abort("All columns in {.arg G} and {.arg W} must be in {.arg p_rwg}.")
        if (ncol(p_rwg) != nlevels(R_vec) + ncol(GW))
            cli_abort("Number of racial categories in {.arg p_rwg} and {.arg R} must match.")
        if (nrow(anti_join(GW, p_rwg, by=colnames(GW))) > 0)
            cli_abort("Some {.arg G}/{.arg W} combinations are missing from {.arg p_rwg}.")
    }
    p_g = prop.table(table(G_vec))
    p_wgr <<- as.matrix(p_rwg[, -1])
    for (i in seq_along(p_r)) {
        p_wgr[, i] = p_wgr[, i] * p_g
        p_wgr[, i] = p_wgr[, i] / sum(p_wgr[, i])
    }

    # TODO remove: debug only
    pp <<- list(
        sr = `rownames<-`(p_sr, levels(S_vec)),
        rs = p_rs,
        wgr = `rownames<-`(p_wgr, levels(G_vec)),
        rwg = p_rwg,
        s = `names<-`(p_s, levels(S_vec)),
        g = `names<-`(p_g, levels(G_vec))
    )

    # baseline
    m_baseline = matrix(nrow=length(X_vec), ncol=length(p_r))
    for (i in seq_along(p_r)) {
        m_baseline[, i] = p_sr[S_vec, i] * p_wgr[G_vec, i] * p_r[i]
    }
    m_baseline = m_baseline / rowSums(m_baseline)
    colnames(m_baseline) = names(p_r)

    m_baseline_surn = matrix(nrow=length(X_vec), ncol=length(p_r))
    for (i in seq_along(p_r)) {
        m_baseline_surn[, i] = p_sr[S_vec, i] * p_r[i]
    }
    m_baseline_surn = m_baseline_surn / rowSums(m_baseline_surn)
    colnames(m_baseline_surn) = names(p_r)

    # inversion
    rlang::check_installed(c("MASS", "glmnet"), "generalized inverse")
    p_rs = p_sr %*% diag(p_r)
    p_rs_inv = t(MASS::ginv(p_rs))
    p_xs = prop.table(table(X_vec, S_vec))

    p_x = prop.table(table(X_vec))
    P_xr_low = outer(p_x, p_r, \(x, y) pmax(0, x+y-1))
    P_xr_high = outer(p_x, p_r, pmin)

    tidy_pred = function(est) {
        est[est > P_xr_high] = P_xr_high[est > P_xr_high]
        est[est < P_xr_low] = P_xr_low[est < P_xr_low]
        rownames(est) = levels(X_vec)
        colnames(est) = names(p_r)
        est
    }

    P_xr_est1 = tidy_pred(p_xs %*% p_rs_inv)
    P_xr_est2 = tidy_pred(rbind(
        as.numeric(glmnet::glmnet(p_rs, t(p_xs[1, ]), lambda=0, lower.limits=0, intercept=FALSE)$beta),
        as.numeric(glmnet::glmnet(p_rs, t(p_xs[2, ]), lambda=0, lower.limits=0, intercept=FALSE)$beta),
        as.numeric(glmnet::glmnet(p_rs, t(p_xs[3, ]), lambda=0, lower.limits=0, intercept=FALSE)$beta),
        as.numeric(glmnet::glmnet(p_rs, t(p_xs[4, ]), lambda=0, lower.limits=0, intercept=FALSE)$beta)
    ) %*% diag(p_r))

    out = list(baseline=m_baseline,
         lsq=P_xr_est1, nnls=P_xr_est2,
         base_surn=m_baseline_surn)

    if (gibbs) {
        out$gibbs = gibbs_race(iter, thin,
                               X_vec, S_vec, G_vec,
                               log(p_sr), log(p_wgr), log(p_r),
                               matrix(alpha, length(p_r), nlevels(X_vec)),
                               verbose)
    }


    if (stan) {
        #tally_df = tibble(X=X_vec, S=S_vec, GZ=G_vec) %>%
        #    count(S, GZ)
        #print(tally_df)
        #print(nrow(data))
        #stop("temp")

        stan_data = list(
            N = length(X_vec),
            n_x = nlevels(X_vec),
            n_r = length(p_r),
            n_gz = nlevels(G_vec),
            n_s = nlevels(S_vec),

            X = as.integer(X_vec),
            S = as.integer(S_vec),
            GZ = as.integer(G_vec),

            p_sr = p_sr,
            p_gzr = p_wgr,
            p_r = p_r,
            p_gz = prop.table(table(G_vec)),

            n_prior_obs = alpha[1]
        )
        #out$stan = rstan::sampling(stanmodels$model_distr, stan_data,
        #                           chains=1, pars="p_xr", iter=500, warmup=400)
        out$stan = rstan::vb(stanmodels$model_distr, stan_data,
                             pars="p_xr", algorithm="meanfield",
                             eta=1, adapt_engaged=FALSE,
                             grad_samples=2, eval_elbo=50,
                             importance_resampling=FALSE)
        #out$stan = rstan::optimizing(stanmodels$model_distr, stan_data,
        #                             verbose=verbose,
        #                             #draws=iter, importance_resampling=TRUE,
        #                             tol_grad=1e-6, tol_param=1e-5, init_alpha=1)
    }

    out
}

# Helper function to make an R|S table
census_surname_table = function(S, S_name, p_r) {
    if (length(p_r) != 5)
        cli_abort(c("Number of racial categories doesn't match the Census table.",
                    "i"="Categories should be White, Black, Hispanic, Asian, and other."))
    require("wru")
    x = data.frame(surname = unique(S)) %>%
        wru::merge_surnames(impute.missing=FALSE) %>%
        suppressWarnings() %>%
        na.omit

    #x$surname = as.character(x$surname)
    x = rbind(x[, -2], list(surname="<generic>", p_whi=p_r[1], p_bla=p_r[2],
                            p_his=p_r[3], p_asi=p_r[4], p_oth=p_r[5]))
    #x = x[match(levels(S), x$surname), ]
    colnames(x) = c(S_name, "white", "black", "hisp", "asian", "other")
    as_tibble(x)
}

# Helper function to make an R|G table
census_zip_table = function(G, G_name, p_r) {
    if (!rlang::is_installed("zipWRUext2"))
        cli_abort(c("{.pkg zipWRUext2} must be installed to use ZIP code information automatically.",
                    ">"=' {.code devtools::install_github("https://github.com/jcuriel-unc/zipWRUext",
                    subdir="zipWRUext2")} to install.'))
    if (length(p_r) != 5)
        cli_abort(c("Number of racial categories doesn't match the Census data.",
                    "i"="Categories should be White, Black, Hispanic, Asian, and other."))

    x = zipWRUext2::zip_all_census2[, c(2, 9:13)]
    x = rbind(x, list(zcta5="<none>", r_whi=p_r[1], r_bla=p_r[2],
                      r_his=p_r[3], r_asi=p_r[4], r_oth=p_r[5]))
    match_idx = match(levels(G), x$zcta5)
    match_idx[is.na(match_idx)] = match("<none>", x$zcta5)
    x = x[match_idx, ]
    zip_sums = rowSums(x[, -1])
    for (i in seq_along(p_r)) {
        x[, 1+i] = x[, 1+i] / zip_sums
    }
    colnames(x) = c(G_name, "white", "black", "hisp", "asian", "other")
    as_tibble(x)
}


#' @export
rake = function(est, row=rowSums(est), col=colSums(est)) {
    for (i in 1:5) {
        est = est %*% diag(col / colSums(est))
        est = diag(row / rowSums(est)) %*% est
    }
    est
}
