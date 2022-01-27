
#' Infer Race Given Name, Location, Outcome, and Covariates
#'
#' Collapsed Gibbs sampler for individual race given name, location, outcome,
#' and covariates.
#'
#' @param X the column containing the outcome
#' @param S the column containing surnames
#' @param G the column containing locations
#' @param Z the column(s), if any, containing other covariates. Use `c()` to provide multiple columns
#' @param data the data
#' @param p_rs a data frame, containing a column matching `S` and columns for
#'   each value of `R` giving the conditional probabilities of R given S.
#'   Defaults to a table from the U.S. Census Bureau from 2010.
#' @param p_rgz a data frame, containing columns matching `G`and `W`, and columns for
#'   each value of `R` giving the conditional probabilities of R given W and G.
#'   Defaults to a table from the 2020 decennial census.
#' @param p_r a vector containing the marginal probabilities for each value of
#'   `R`. Defaults to the demographics of North Carolina.
#' @param regularize if `TRUE`, regularize the Census-calculated `R|S` and `R|G,Z` tables
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
model_race = function(X, S, G, Z=NULL, data=NULL, p_rs=NULL, p_rgz=NULL,
                      p_r=c(white=0.626, black=0.222, hisp=0.098, asian=0.032, other=0.022),
                      regularize=TRUE, alpha=7,
                      methods=c("bis", "bisg", "gibbs", "lsq", "nonparam", "additive"),
                      stan_method=c("vb", "opt", "hmc"),
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
    Z_df = data[, eval_select(enquo(Z), data)]

    check_vec = function(x) (is.character(x) | is.factor(x)) && !any(is.na(x))

    if (!check_vec(X_vec)) cli_abort("{.arg X} must be a character or factor with no missing values.")
    if (!check_vec(S_vec)) cli_abort("{.arg S} must be a character or factor with no missing values.")
    if (!is.character(G_vec) && !is.factor(G_vec))
        cli_abort("{.arg G} must be a character or factor vector.")
    if (any(is.na(G_vec))) {
        #cli_inform("Missing values found in {.arg G}")
        G_vec = coalesce(G_vec, "<none>")
    }
    if (!all(vapply(Z_df, class, character(1)) %in% c("character", "factor"))) {
        cli_abort("{.arg Z} must contain only character or factor columns.")
    }
    if (any(is.na(Z_df))) cli_abort("Missing values found in {.arg Z}")

    X_vec = as.factor(X_vec)
    S_vec = as.factor(S_vec)
    G_vec = as.factor(G_vec)
    GZ = cbind(G_vec, Z_df)
    GZ_vec = as.factor(vctrs::vec_duplicate_id(GZ))
    GZ_levels = vapply(GZ, nlevels, integer(1))
    n_gz_col = sum(GZ_levels)
    GZ_var = inverse.rle(list(lengths=GZ_levels, values=1:ncol(GZ)))
    GZ_mat = do.call(cbind, lapply(GZ, function(x) {
        out = matrix(0, nrow=length(x), ncol=nlevels(x))
        out[cbind(seq_along(x), as.integer(x))] = 1
        out
    }))

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

    if (missing(p_rgz) & !missing(Z)) { # TODO remove this
        names(GZ)[1] = "zip"
        p_rgz = census_zip_age_sex_table(GZ, GZ_vec, p_r, regularize)
    } else if (missing(p_rgz)) {
        p_rgz = census_zip_table(G_vec, as_name(enquo(G)), p_r, regularize)
    } else {
        if (!is.data.frame(p_rs)) cli_abort("{.arg p_rgz} must be a data frame.")
        if (!all(colnames(GZ) %in% colnames(p_rgz)))
            cli_abort("All columns in {.arg G} and {.arg Z} must be in {.arg p_rgz}.")
        if (ncol(p_rgz) != nlevels(R_vec) + ncol(GZ))
            cli_abort("Number of racial categories in {.arg p_rgz} and {.arg R} must match.")
        if (nrow(anti_join(GW, p_rgz, by=colnames(GZ))) > 0)
            cli_abort("Some {.arg G}/{.arg Z} combinations are missing from {.arg p_rgz}.")
    }
    if (missing(Z)) {
        p_gz = prop.table(table(G_vec))
    } else {
        p_gz = prop.table(table(GZ_vec))
    }
    p_gzr = as.matrix(select(p_rgz, white:other))
    for (i in seq_along(p_r)) {
        p_gzr[, i] = p_gzr[, i] * p_gz
        p_gzr[, i] = p_gzr[, i] / sum(p_gzr[, i])
    }

    out = list()
    methods = match.arg(methods, several.ok=TRUE)

    if ("bis" %in% methods) {
        out$bis = est_bisg(X_vec, S_vec, G_vec, p_sr, p_gzr, p_r, geo=FALSE)
    }
    pr_base = est_bisg(X_vec, S_vec, GZ_vec, p_sr, p_gzr, p_r, geo=TRUE)
    if ("bisg" %in% methods) {
        out$bisg = pr_base
    }

    if ("lsq" %in% methods) {
        out = c(out, est_leastsq(X_vec, S_vec, p_sr, p_r))
    }

    if ("gibbs" %in% methods) {
        M_sr = census_surname_table(S_vec, as_name(enquo(S)), p_r, FALSE, count=TRUE)[, -1] %>%
            as.matrix()
        if (!missing(Z)) {
            N_gzr = census_zip_age_sex_table(GZ, GZ_vec, p_r, FALSE, count=TRUE) %>%
                select(white:other) %>%
                as.matrix()
        } else {
            N_gzr = census_zip_table(G_vec, as_name(enquo(G)), p_r, FALSE, count=TRUE)[, -1] %>%
                as.matrix()
        }
        alpha_sr = matrix(1/nrow(M_sr), nrow(M_sr), ncol(M_sr))
        beta_gzr = matrix(rep(p_r, each=nrow(N_gzr)), nrow=nrow(N_gzr))

        out$gibbs = gibbs_race(iter, min(100, floor(iter/2)), X_vec, S_vec, GZ_vec,
                               M_sr, N_gzr, alpha_sr, beta_gzr, verbose)
    }

    # TODO remove
    d <<- list(X=X_vec, GZ=GZ_vec, GZ_mat=GZ_mat, GZ_var=GZ_var, pr_base=pr_base)

    if ("nonparam" %in% methods) {
        out$nonparam = est_nonparam(X_vec, GZ_vec, pr_base, alpha[1],
                                    match.arg(stan_method), iter, verbose)
    }
    gc()
    if ("additive" %in% methods) {
        out$additive = est_additive(X_vec, GZ_mat, GZ_var, pr_base,
                                    match.arg(stan_method), iter, verbose)
    }
    gc()

    out
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
