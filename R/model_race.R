
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
#'
#' @return TBD
#' @export
model_race = function(X, S, G, W=NULL, data=NULL, p_rs=NULL, p_rwg=NULL,
                      p_r=c(white=0.626, black=0.222, hisp=0.098, asian=0.032, other=0.022),
                      alpha=1,
                      iter=100, thin=1) {
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

    ## Parse anc check input probabilities ----------------
    if (missing(p_rs)) {
        S_vec = as.character(S_vec)
        S_vec[!S_vec %in% wru::surnames2010$surname] = "<generic>"
        S_vec = as.factor(S_vec)
        p_rs = census_surname_table(S_vec, as_name(enquo(S)), p_r)
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
    p_wgr = as.matrix(p_rwg[, -1])
    for (i in seq_along(p_r)) {
        p_wgr[, i] = p_wgr[, i] * p_g
        p_wgr[, i] = p_wgr[, i] / sum(p_wgr[, i])
    }

    gibbs_race(iter, thin,
               as.integer(X_vec), as.integer(S_vec), as.integer(G_vec),
               log(p_sr), log(p_wgr), log(p_r),
               matrix(alpha, nlevels(X_vec), length(p_r)))
}

if (F) {
model_race(party, last_name, zip, data=voters)
}

# Helper function to make an R|S table
census_surname_table = function(S, S_name, p_r) {
    x = wru::surnames2010
    if (ncol(x) - 1 != length(p_r))
        cli_abort(c("Number of racial categories doesn't match the Census table.",
                    "i"="Categories should be White, Black, Hispanic, Asian, and other."))

    x$surname = as.character(x$surname)
    x = rbind(x, list(surname="<generic>", p_whi=p_r[1], p_bla=p_r[2],
                      p_his=p_r[3], p_asi=p_r[4], p_oth=p_r[5]))
    x = x[match(levels(S), x$surname), ]
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
    x = x[match(levels(G), x$zcta5), ]
    colnames(x) = c(G_name, "white", "black", "hisp", "asian", "other")
    as_tibble(x)
}

