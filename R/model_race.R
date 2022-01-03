
#' Infer Race Given Name, Location, Outcome, and Covariates
#'
#' @param X the column containing the outcome
#' @param R the column containing race labels
#' @param S the column containing surnames
#' @param G the column containing locations
#' @param W the column(s), if any, containing other covariates. Use `c()` to provide multiple columns
#' @param data the data
#' @param p_rs a data frame, containing a column matching `S` and columns for
#'   each value of `R` giving the conditional probabilities of R given S.
#'   Defaults to a table from the U.S. Census Bureau from 2010.
#' @param p_wgr a data frame, containing columns matching `G`and `W`, and columns for
#'   each value of `R` giving the conditional probabilities of R given W and G.
#'   Defaults to a table from the 2020 decennial census.
#' @param p_r a vector containing the marginal probabilities for each value of
#'   `R`. Defaults to the demographics of North Carolina.
#'
#' @return TBD
#' @export
model_race = function(X, R, S, G, W=NULL, data=NULL, p_rs=NULL, p_wgr=NULL,
                      p_r=c(white=0.626, black=0.222, hisp=0.098, asian=0.032, other=0.022)) {
    ## Parse and check input data frame ----------------
    if (missing(data)) cli_abort("{.arg data} must be provided.")
    X_vec = eval_tidy(enquo(X), data)
    R_vec = eval_tidy(enquo(R), data)
    S_vec = eval_tidy(enquo(S), data)
    G_vec = eval_tidy(enquo(G), data)
    W_df = data[, eval_select(enquo(W), data)]

    check_vec = function(x) (is.character(x) | is.factor(x)) && !any(is.na(x))

    if (!check_vec(X_vec)) cli_abort("{.arg X} must be a character or factor with no missing values.")
    if (!check_vec(R_vec)) cli_abort("{.arg R} must be a character or factor with no missing values.")
    if (!check_vec(S_vec)) cli_abort("{.arg S} must be a character or factor with no missing values.")
    if (!is.character(G_vec) && !is.factor(G_vec))
        cli_abort("{.arg G} must be a character or factor vector.")
    if (any(is.na(G_vec))) {
        cli_inform("Missing values found in {.arg G}")
        G_vec = coalesce(G_vec, "<none>")
    }
    if (!all(vapply(W_df, class, character(1)) %in% c("character", "vector"))) {
        cli_abort("{.arg W} must contain only character or factor columns.")
    }
    if (any(is.na(W_df))) cli_abort("Missing values found in {.arg W}")

    X_vec = as.factor(X_vec)
    R_vec = as.factor(R_vec)
    S_vec = as.factor(S_vec)
    GW = cbind(G_vec, W_df)

    ## Parse anc check input probabilities ----------------
    if (missing(p_rs)) {
        p_rs = census_surname_table(R_vec, as_name(S))
    } else {
        if (!is.data.frame(p_rs)) cli_abort("{.arg p_rs} must be a data frame.")
        if (ncol(p_rs) != nlevels(R_vec) + 1L)
            cli_abort("Number of racial categories in {.arg p_rs} and {.arg R} must match.")
        if (!all(S_vec %in% p_rs[, 1]))
            cli_abort("Some names are missing from {.arg p_rs}.")
    }
    # TODO convert R|S to S|R

    if (missing(p_wgr)) {
        p_rs = census_zip_table(R_vec, G_vec, p_r)
    } else {
        if (!is.data.frame(p_rs)) cli_abort("{.arg p_wgr} must be a data frame.")
        if (!all(colnames(GW) %in% colnames(p_wgr)))
            cli_abort("All columns in {.arg G} and {.arg W} must be in {.arg p_wgr}.")
        if (ncol(p_wgr) != nlevels(R_vec) + ncol(GW))
            cli_abort("Number of racial categories in {.arg p_wgr} and {.arg R} must match.")
        if (nrow(anti_join(GW, p_wgr, by=colnames(GW))) > 0)
            cli_abort("Some {.arg G}/{.arg W} combinations are missing from {.arg p_wgr}.")
    }
    # TODO convert R|W,G to W,G|R


    invisible(TRUE)
}

if (F) {
model_race(party, race, last_name, zip, data=voters)
}

# Helper function to make an R|S table
census_surname_table = function(R, S_name) {
    x = wru::surnames2010
    if (ncol(x) - 1 != nlevels(R))
        cli_abort(c("Number of racial categories doesn't match the Census table.",
                    "i"="Categories should be White, Black, Hispanic, Asian, and other."))
    colnames(x) = c(S_name, "white", "black", "hisp", "asian", "other")
    x
}

# Helper function to make an R|G table
census_zip_table = function(R, G, p_r) {
    rlang::check_installed("zipWRUext2", "for ZIP/race data.")
    x = zipWRUext2::zip_all_census2
    if (nlevels(R) != 5)
        cli_abort(c("Number of racial categories doesn't match the Census data.",
                    "i"="Categories should be White, Black, Hispanic, Asian, and other."))
    colnames(x)[1] = c(as_name(enquo(S)),
                       "white", "black", "hisp", "asian", "other")
    x
}

