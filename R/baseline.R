
#' Infer Race Given Name, Location, Outcome, and Covariates
#'
#' BISG+ race probability calculations.
#'
#' @param X the column containing the outcome
#' @param S the column containing surnames
#' @param G the column containing ZIP codes
#' @param Z the column(s), if any, containing other covariates. Use `c()` to provide multiple columns
#' @param data the data
#' @param p_rs a data frame, containing a column matching `S` and columns for
#'   each value of `R` giving the conditional probabilities of R given S.
#'   Defaults to a table from the U.S. Census Bureau from 2010.
#' @param p_rgz a data frame, containing columns matching `G`and `Z`, and columns for
#'   each value of `R` giving the conditional probabilities of R given Z and G.
#'   Defaults to a table from the 2020 decennial census.
#' @param p_r a vector containing the marginal probabilities for each value of
#'   `R`. Defaults to the demographics of the US.
#' @param iterate how much to iteratively refine the race prior from the data. Set to 0 to disable.
#' @param regularize if `TRUE`, regularize the Census-calculated `R|S` and `R|G,Z` tables
#'
#' @return a tibble, with rows matching `data` and columns for the race probabilities
#' @export
predict_race_sgz = function(S, G, Z=NULL, data=NULL, p_rs=NULL, p_rgz=NULL,
                            p_r=c(white=0.615, black=0.123, hisp=0.176, asian=0.053, other=0.034),
                            iterate=10L, regularize=TRUE) {
    ## Parse and check input data frame ----------------
    if (missing(data)) cli_abort("{.arg data} must be provided.")
    S_vec = eval_tidy(enquo(S), data)
    G_vec = eval_tidy(enquo(G), data)
    Z_df = data[, eval_select(enquo(Z), data)]

    if (!check_vec(S_vec)) cli_abort("{.arg S} must be a character or factor with no missing values.")
    if (!is.character(G_vec) && !is.factor(G_vec))
        cli_abort("{.arg G} must be a character or factor vector.")
    if (!all(vapply(Z_df, class, character(1)) %in% c("character", "factor"))) {
        cli_abort("{.arg Z} must contain only character or factor columns.")
    }
    if (any(is.na(Z_df))) cli_abort("Missing values found in {.arg Z}")

    G_vec = as.factor(coalesce(G_vec, "<none>"))
    GZ = cbind(G_vec, Z_df)
    GZ_vec = as.factor(vctrs::vec_duplicate_id(GZ))

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

    if (missing(p_rgz) && !missing(Z) && ncol(Z_df) == 2) { # TODO remove this
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
        if (nrow(anti_join(GZ, p_rgz, by=colnames(GZ))) > 0)
            cli_abort("Some {.arg G}/{.arg Z} combinations are missing from {.arg p_rgz}.")
    }
    if (missing(Z)) {
        GZ_vec = G_vec
    }
    p_gz = prop.table(table(GZ_vec))

    # flip which margin sums to 1
    p_gzr = as.matrix(select(p_rgz, white:other))
    for (i in seq_along(p_r)) {
        p_gzr[, i] = p_gzr[, i] * p_gz
        p_gzr[, i] = p_gzr[, i] / sum(p_gzr[, i])
    }

    m_bisg = est_bisg(S_vec, GZ_vec, p_sr, p_gzr, p_r)
    if (iterate > 0) {
        for (i in 1:iterate) {
            p_r[1:5] = colMeans(m_bisg)
            m_bisg = est_bisg(S_vec, GZ_vec, p_sr, p_gzr, p_r)
        }
    }

    list(pr=as_tibble(m_bisg), GZ=GZ_vec, S=p_sr[S_vec, ])
}


est_bisg = function(S, G, p_sr, p_gzr, p_r, geo=TRUE) {
    if (!geo) p_gzr = p_gzr*0 + 1

    m_bisg = matrix(nrow=length(S), ncol=length(p_r))
    for (i in seq_along(p_r)) {
        m_bisg[, i] = p_sr[S, i] * p_gzr[G, i] * p_r[i]
    }
    m_bisg = m_bisg / rowSums(m_bisg)
    colnames(m_bisg) = paste0("pr_", names(p_r))
    m_bisg
}
