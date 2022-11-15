
#' Infer Race Given Name, Location, Outcome, and Covariates
#'
#' BISG(Z) race probability calculations.
#'
#' @param S the column containing surnames
#' @param G the column containing ZIP codes
#' @param Z the column(s), if any, containing other covariates. Use `c()` to provide multiple columns
#' @param data the data
#' @param p_rs a data frame, containing a column matching `S` and columns for
#'   each value of `R` giving the conditional probabilities of R given S.
#'   Defaults to a table from the U.S. Census Bureau from 2010.
#' @param p_rgz a data frame, containing columns matching `G`and `Z`, and columns for
#'   each value of `R` giving the conditional probabilities of R given Z and G.
#'   Defaults to a ZIP code table from the 2010 decennial census.
#' @param p_r a vector containing the marginal probabilities for each value of
#'   `R`. Defaults to the demographics of the US if `est_r_gz` is `FALSE`.
#' @param est_r_gz If `TRUE` estimate race prior `p_r` from distribution of `G` and `Z` in sample
#' @param iterate how much to iteratively refine the race prior from the data. Set to 0 to disable.
#' @param return_gzr whether to return an estimated p(G,Z|R) matrix as the
#'   `p_gzr` attribute and combined G,Z vectors as the `gz` attribute.
#'   Set to `FALSE` to save space if these values are not needed.
#'   They are required to use [calc_joint_bisgz_ols]
#'
#' @return a tibble, with rows matching `data` and columns for the race probabilities
#' @export
predict_race_sgz = function(S, G, Z=NULL, data=NULL, p_rs=NULL, p_rgz=NULL,
                            p_r=c(white=0.630, black=0.121, hisp=0.173,
                                  asian=0.0478, aian=0.0072, other=0.0210),
                            est_r_gz=TRUE, iterate=8L * (!est_r_gz & nrow(data) > 100),
                            return_gzr=TRUE) {
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

    G_name = as_name(enquo(G))
    G_vec = as.factor(coalesce(G_vec, "<none>"))
    GZ = dplyr::bind_cols("{G_name}":=G_vec, Z_df)
    GZ_vec = as.factor(vctrs::vec_duplicate_id(GZ))

    ## Parse and check input probabilities ----------------
    if (missing(p_rs)) {
        S_vec = as.character(S_vec)
        p_rs = census_surname_table(S_vec, "last_name", p_r, flip=TRUE)
        S_vec[!S_vec %in% p_rs[[1]]] = "<generic>"
        S_vec = factor(S_vec, levels=p_rs[[1]])
        p_sr = as.matrix(p_rs[, -1])
    } else {
        if (!is.data.frame(p_rs)) cli_abort("{.arg p_rs} must be a data frame.")
        if (!all(S_vec %in% p_rs[, 1]))
            cli_abort("Some names are missing from {.arg p_rs}.")

        p_s = prop.table(table(S_vec))
        p_sr = as.matrix(p_rs[, -1])
        for (i in seq_along(p_r)) {
            p_sr[, i] = p_sr[, i] * p_s
            p_sr[, i] = p_sr[, i] / sum(p_sr[, i])
        }
    }

    if (missing(p_rgz) && !missing(Z) && ncol(Z_df) == 2) { # TODO remove this (NC-specific)
        names(GZ)[1] = "zip"
        #p_rgz = census_zip_age_sex_table(GZ, GZ_vec, p_r, regularize)
    } else if (missing(p_rgz)) {
        p_rgz = census_zip_table(G_vec, G_name, p_r, counts=TRUE)
    } else {
        if (!is.data.frame(p_rs)) cli_abort("{.arg p_rgz} must be a data frame.")
        if (!all(colnames(GZ) %in% colnames(p_rgz)))
            cli_abort("All columns in {.arg G} and {.arg Z} must be in {.arg p_rgz}.")
        if (ncol(p_rgz) != 6 + ncol(GZ))
            cli_abort("Number of racial categories in {.arg p_rgz} and {.arg R} must match.")
        # if (nrow(anti_join(GZ, p_rgz, by=colnames(GZ))) > 0)
        #     cli_abort("Some {.arg G}/{.arg Z} combinations are missing from {.arg p_rgz}.")
    }
    if (missing(Z)) {
        GZ_vec = G_vec
    }

    # match geographies
    if (!"<none>" %in% p_rgz[[G_name]]) {
        p_rgz = rbind(p_rgz, rlang::list2("{G_name}":="<none>",
                      white=p_r[1], black=p_r[2], hisp=p_r[3],
                      asian=p_r[4], aian=p_r[5], other=p_r[6]))
    }
    p_rgz = mutate(p_rgz, across(.data$white:.data$other, ~ coalesce(., p_r[cur_column()])))
    match_idx = match(levels(G_vec), p_rgz[[G_name]])
    match_idx[is.na(match_idx)] = match("<none>", p_rgz[[G_name]])
    p_rgz = p_rgz[match_idx, ]

    # flip which margin sums to 1
    p_gz = prop.table(table(GZ_vec))
    p_gzr = as.matrix(select(p_rgz, .data$white:.data$other))
    for (i in seq_along(p_r)) {
        p_gzr[, i] = p_gzr[, i] * p_gz
        p_gzr[, i] = p_gzr[, i] / sum(p_gzr[, i])
    }

    if (isTRUE(est_r_gz)) {
        p_r = colSums(p_rgz[, -1])
        p_r = p_r / sum(p_r)
    }

    m_bisg = est_bisg(S_vec, GZ_vec, p_sr, p_gzr, p_r)
    if (iterate > 0) {
        for (i in 1:iterate) {
            p_r[1:6] = colMeans(m_bisg)
            m_bisg = est_bisg(S_vec, GZ_vec, p_sr, p_gzr, p_r)
        }
    }

    out = as_tibble(m_bisg)
    if (isTRUE(return_gzr)) {
        rownames(p_gzr) = p_rgz[[G_name]]
        attr(out, "p_gzr") = p_gzr
        attr(out, "gz") = GZ_vec
    }
    out
}


est_bisg = function(S, GX, p_sr, p_gxr, p_r, geo=TRUE) {
    if (!geo) p_gxr = p_gxr*0 + 1

    m_bisg = calc_bayes_bisg(S, GX, p_sr, p_gxr, p_r);
    colnames(m_bisg) = paste0("pr_", names(p_r))
    m_bisg
}
