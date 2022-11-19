#' Bayesian Improved Surname Geocoding (BISG)
#'
#' Calculates individual probabilities of belonging to racial groups given last
#' name, location, and other covariates (optional). The standard function
#' `bisg()` treats the input tables as fixed. An alternative function
#' `bisg_me()`, assumes that the input tables are subject to measurement error,
#' and uses a Gibbs sampler to impute the individual race probabilities, using
#' the model of Imai et al. (2022).
#'
#' @param formula A formula specifying the BISG model. Must include the special
#'   term `nm()` to identify the surname variable. Certain geographic variables
#'   can be identified similarly: `zip()` for ZIP codes, and `county()` for
#'   counties. If no other predictor variables are provided, then `bisg()` will
#'   automatically be able to build a table of census data to use in inference.
#'   If other predictor variables are included, or if other geographic
#'   indicators are used, then the user must specify the `p_rgx` argument below.
#'   The left-hand side of the formula is ignored.
#'   See the examples section below for sample formulas.
#' @param data The data frame containing the variables in `formula`.
#' @param p_r The prior distribution of race in the sample, as a numeric vector.
#'   Defaults to U.S. demographics as provided by [p_r_natl()].
#'   Can also set `p_r="est"` or `"estimate"` to estimate this from the
#'   geographic distribution.  Since the prior distribution on race strongly
#'   affects the calibration of the BISG probabilities and thus the accuracy of
#'   downstream estimates, users are encouraged to think carefully about an
#'   appropriate value for `p_r`.  If no prior information on the racial makeup
#'   of the sample is available, and yet the sample is very different from the
#'   overall U.S. population, then `p_r="estimate"` will likely produce superior
#'   results.
#' @param p_rgx The distribution of race given location (G) and other covariates
#'   (X) specified in `formula`. Should be provided as a data frame, with columns
#'   matching the predictors in `formula`, and additional columns for each
#'   racial group containing the conditional probability for that racial group
#'   given the predictors. For example, if Census tracts are the only predictors,
#'   `p_rgx` should be a data frame with a tract column and columns `white`,
#'   `black`, etc. containing the racial distribution of each tract.
#'   If `formula` contains only labeled terms (like `zip()`), then by default
#'   `p_rgx` will be constructed automatically from the most recent Census data.
#'   This table will be normalized by row, so it can be provided as population
#'   counts as well. Counts are required for `bisg_me()`.
#'   TODO LINK TO CENSUS FUNCTIONS
#' @param p_rs The distribution of race given last name. As with `p_rgx`, should
#'   be provided as a data frame, with a column of names and additional columns
#'   for each racial group. Users should not have to specify this argument in
#'   most cases, as the table will be built from published Census surname tables
#'   automatically. Counts are required for `bisg_me()`.
#'
#' @return An object of class `bisg`, which is just a data frame with some
#'   additional attributes. The data frame has rows matching the input data and
#'   columns for the race probabilities.
#'
#' @examples
#' data(pseudo_vf)
#' bisg(~ nm(last_name), data=pseudo_vf)
#' bisg(~ nm(last_name) + zip(zip), data=pseudo_vf)
#'
#' @references
#' Elliott, M. N., Fremont, A., Morrison, P. A., Pantoja, P., and Lurie, N.
#' (2008). A new method for estimating race/ethnicity and associated disparities
#' where administrative records lack self-reported race/ethnicity. *Health
#' services research*, 43(5p1):1722–1736.
#'
#' Fiscella, K. and Fremont, A. M. (2006). Use of geocoding and surname analysis
#' to estimate race and ethnicity. *Health services research*,
#' 41(4p1):1482–1500.
#'
#' Imai, K., Olivella, S., and Rosenman, E. T. (2022). Addressing census data
#' problems in race imputation via fully Bayesian improved surname geocoding and
#' name supplements. *arXiv preprint arXiv:2205.06129*.
#'
#' @describeIn bisg The standard BISG model.
#' @export
bisg <- function(formula, data=NULL, p_r=p_r_natl(), p_rgx=NULL, p_rs=NULL) {
    vars = parse_bisg_form(formula, data)

    l_name = make_name_tbl_vec(vars, p_r, p_rs, FALSE)
    l_gx = make_gx_tbl_vec(vars, p_r, p_rgx)

    m_bisg = est_bisg(l_name$S, l_gx$GX, l_name$p_sr, l_gx$p_gxr, l_gx$p_r)

    out <- as_tibble(m_bisg)
    class(out) = c("bisg", class(out))
    attr(out, "S_name") = vars$S_name
    attr(out, "p_r") = l_gx$p_r
    attr(out, "method") = "std"

    out
}

#' @param iter How many sampling iterations in the Gibbs sampler
#' @param warmup How many burn-in iterations in the Gibbs sampler
#'
#' @examples
#' data(pseudo_vf)
#' bisg_me(~ nm(last_name) + zip(zip), data=pseudo_vf)

#' @describeIn bisg The measurement error BISG model.
#' @export
bisg_me <- function(formula, data=NULL, p_r=p_r_natl(), p_rgx=NULL, p_rs=NULL,
                    iter=1000, warmup=100) {
    vars = parse_bisg_form(formula, data)

    l_name = make_name_tbl_vec(vars, p_r, p_rs, TRUE)
    l_gx = make_gx_tbl_vec(vars, p_r, p_rgx)

    n_gx = nrow(l_gx$p_rgx)
    n_s = nrow(l_name$p_sr)
    alpha_gzr = matrix(rep(l_gx$p_r, n_gx), nrow=n_gx, ncol=6, byrow=T)
    beta_sr = matrix(rep(l_gx$p_r, n_s), nrow=n_s, ncol=6, byrow=T)

    m_bisg = gibbs_me(iter+warmup, warmup, l_name$S, l_gx$GX,
                      l_name$p_sr, l_gx$p_rgx,
                      alpha_gzr, beta_sr, verbosity=3L)
    colnames(m_bisg) = paste0("pr_", names(p_r))

    out <- as_tibble(m_bisg)
    class(out) = c("bisg", class(out))
    attr(out, "S_name") = vars$S_name
    attr(out, "p_r") = l_gx$p_r
    attr(out, "method") = "me"

    out
}

# Parse formula and return a list with name vector, name vector name,
# type of geography used (or "none"), and data frame of other covariates
parse_bisg_form <- function(formula, data=NULL) {
    formula = update.formula(formula, NULL ~ .)
    f_terms = terms(formula, specials=c("nm", "zip", "county"), data=data)
    d_model = get_all_vars(f_terms, data=data)
    if (ncol(d_model) != length(attr(f_terms, "term.labels"))) {
        cli_abort("Duplicated variables found in {.arg formula}.", call=parent.frame())
    }

    # check name
    nm_loc = attr(f_terms, "specials")$nm
    if (is.null(nm_loc)) {
        cli_abort("{.arg formula} must have a {.fn nm} term
                  specifying the surname variable.", call=parent.frame())
    }
    S = d_model[[nm_loc]]
    if (!check_vec(S)) {
        cli_abort("The names vector provided in {.fn nm} must be
                   a character or factor with no missing values.",
                  call=parent.frame())
    }

    # check ZIP, county, etc.
    geo_locs = vapply(attr(f_terms, "specials")[-1], function(x) {
        if (is.null(x)) NA_integer_ else x
    }, 1L)
    n_loc = sum(!is.na(geo_locs))
    geo_type = "none"
    if (n_loc > 1) {
        cli_abort("{.arg formula} can only have one labeled location vector.
                  {.fn {names(which(!is.na(geo_locs)))}} were found.",
                  call=parent.frame())
    } else if (n_loc == 1 && ncol(d_model) > 2) {
        cli_warn("Labeled location vectors in {.arg formula} are only used if
                 there are no other predictors.", call=parent.frame())
    } else if (n_loc == 1) {
        geo_type = names(which(!is.na(geo_locs)))
        geo_loc = geo_locs[geo_type]
        if (geo_type == "county") {
            cli_abort("{.fn county} is not yet supported.", call=parent.frame())
        }

        d_model[[geo_loc]] =  as.factor(coalesce(d_model[[geo_loc]], "<none>"))
    }

    GX = d_model[-nm_loc]
    if (!all(vapply(GX, class, character(1)) %in% c("character", "factor"))) {
        cli_abort("Predictors in {.arg formula} must be character or factor vectors.",
                  call=parent.frame())
    }
    if (ncol(GX) > n_loc && any(is.na(GX))) {
        cli_abort("Missing values found in {.arg formula} predictors",
                  call=parent.frame())
    }

    list(S = S,
         S_name = colnames(d_model)[nm_loc],
         geo_type = geo_type,
         GX = GX)
}

# Prepare name vector and P(S | R) table
make_name_tbl_vec <- function(vars, p_r, p_rs, for_me=FALSE) {
    S = as.character(vars$S)
    if (is.null(p_rs)) {
        S = proc_names(S)
        if (is.character(p_r)) p_r = p_r_natl()

        if (isTRUE(for_me)) {
            p_rs = census_surname_table(S, vars$S_name, p_r, counts=TRUE, flip=FALSE)
        } else {
            p_rs = census_surname_table(S, vars$S_name, p_r, counts=FALSE, flip=TRUE)
        }

        S[!S %in% p_rs[[1]]] = "<generic>" # p_rs[[1]] is the name column
        S = factor(S, levels=p_rs[[1]])
        p_sr = as.matrix(p_rs[, -1])
    } else {
        if (!is.data.frame(p_rs)) {
            cli_abort("{.arg p_rs} must be a data frame.", call=parent.frame())
        }
        if (anyDuplicated(p_rs) > 0) {
            cli_abort("{.arg p_rs} must have unique rows.", call=parent.frame())
        }

        name_col = match(vars$S_name, colnames(p_rs))
        if (is.na(name_col)) {
            cli_abort("Name column {.var {vars$S_name}} not found in {.arg p_rs}.",
                      call=parent.frame())
        }
        if (!all(S %in% p_rs[[name_col]])) {
            cli_abort("Some names are missing from {.arg p_rs}.", call=parent.frame())
        }
        if (any(is.nan(p_rs) | is.na(p_rs) | p_rs < 0 | is.infinite(p_rs))) {
            cli_abort("{.arg p_rs} contains missing, negative,
                      or otherwise invalid values.", call=parent.frame())
        }

        S = factor(S, levels=p_rs[[name_col]])

        p_s = prop.table(table(S))
        p_rs = as.matrix(p_rs[, -name_col])
        p_sr = p_rs / rowSums(p_rs)
        for (i in seq_along(p_r)) {
            p_sr[, i] = p_sr[, i] * p_s
            p_sr[, i] = p_sr[, i] / sum(p_sr[, i])
        }
    }

    list(S = S,
         p_sr = p_sr) # p_sr is actualy p_rs, unormalized, if `for_me=TRUE`
}

# Prepare geo/covariate vector and P(G, X | R) table
make_gx_tbl_vec <- function(vars, p_r, p_rgx) {
    if (ncol(vars$GX) == 0) { # no covariates
        if (is.character(p_r)) {
            cli_abort("{.arg p_r} must be a numeric vector if no
                      other predictors are used.", call=parent.frame())
        }

        N = length(vars$S)
        return(list(GX = rep(1L, N),
                    p_r = p_r,
                    p_rgx = matrix(1e5 * p_r, nrow=1),
                    p_gxr = matrix(rep(1, length(p_r)), nrow=1)))
    }


    if (vars$geo_type == "zip") {
        p_r_tmp = p_r
        if (is.character(p_r)) p_r_tmp = p_r_natl()

        G_name = names(vars$GX)
        p_rgx = census_zip_table(vars$GX[[1]], G_name, p_r_tmp, counts=TRUE)

        # tidy up
        if (!"<none>" %in% p_rgx[[G_name]]) {
            new_row = as.data.frame(as.list(1e5 * p_r_tmp))
            new_row[[G_name]] = "<none>"
            p_rgx = rbind(p_rgx, new_row)
        }
        p_rgx = dplyr::mutate(
            p_rgx, dplyr::across(.data$white:.data$other,
                                 ~ dplyr::coalesce(., p_r_tmp[dplyr::cur_column()]))
        )

        match_idx = match(vars$GX[[1]], p_rgx[[G_name]])
        idx_miss = which(is.na(match_idx))
        if (length(idx_miss) > 0) {
            vars$GX[[1]][idx_miss] = "<none>"
        }
    } else if (vars$geo_type == "county") {
        # TODO IMPLEMENT
    } else { # "none" -> use `p_rgx`
        if (!is.data.frame(p_rgx)) {
            cli_abort("{.arg p_rgx} must be a data frame.", call=parent.frame())
        }
        if (!all(colnames(S$GX) %in% colnames(p_rgx))) {
            cli_abort("All predictor columns must be in {.arg p_rgx}.", call=parent.frame())
        }
        if (!is.character(p_r) && ncol(p_rgx) != length(p_r) + ncol(S$GX)) {
            cli_abort("Number of racial categories in {.arg p_rgx}
                      and {.arg p_r} must match.", call=parent.frame())
        }
        if (anyDuplicated(p_rgx) > 0) {
            cli_abort("{.arg p_rgx} must have unique rows.", call=parent.frame())
        }

        d_miss <- dplyr::anti_join(vars$GX, p_rgx, by=names(vars$GX))
        if (nrow(d_miss) > 0) {
            cli::cli_text("Missing from {.arg p_rs}:")
            print(d_miss)
            cli_abort("Some predictor combinations are missing from {.arg p_rs}.",
                      call=parent.frame())
        }
    }

    # subset to needed rows
    d_match = dplyr::left_join(vars$GX, p_rgx, by=names(vars$GX))
    GX_vec = as.factor(vctrs::vec_duplicate_id(d_match))
    p_rgx = d_match[!duplicated(GX_vec), ]

    # flip which margin sums to 1
    p_gx = prop.table(table(GX_vec))
    p_rgx = as.matrix(dplyr::select(p_rgx, -names(vars$GX)))
    p_gxr = p_rgx / rowSums(p_rgx)
    for (i in seq_len(ncol(p_gxr))) {
        p_gxr[, i] = p_gxr[, i] * p_gx
        p_gxr[, i] = p_gxr[, i] / sum(p_gxr[, i])
    }

    if (is.character(p_r) && (p_r == "est" || p_r == "estimate")) {
        p_r = colSums(p_rgx)
        p_r = p_r / sum(p_r)
    }

    list(GX = GX_vec,
         p_r = p_r,
         p_rgx = p_rgx,
         p_gxr = p_gxr)
}

#' National Racial Demographics
#'
#' Returns the proportion of the U.S. population in six racial groups in a given
#' year. Group definitions necessarily follow those used by the Census Bureau in
#' its surname tables:
#' * `white`: Non-Hispanic White alone
#' * `black`: Non-Hispanic Black alone
#' * `hisp`: Hispanic, any race
#' * `asian`: Non-Hispanic Asian, Native Hawaiian, or Pacific Islander alone
#' * `aian`: Non-Hispanic American Indian/Alaska Native
#' * `other`: Non-Hispanic, two or more races, or other race
#'
#' @param year The year to return demographics for.
#' @param vap If `TRUE`, return statistics for the voting-age population (18+)
#'   rather than the full U.S. population.
#'
#' @returns A named numeric vector of length 6.
#'
#' @examples
#' p_r_natl(year=2010)
#'
#' @export
p_r_natl <- function(year=2010, vap=FALSE) {
    if (isTRUE(vap)) cli_abort("Only {.arg vap = FALSE} is supported for now.")

    year = as.integer(year)
    if (year == 2010L) {
        c(white=0.630, black=0.121, hisp=0.173,
          asian=0.0478, aian=0.0072, other=0.0210)
    } else {
        cli_abort("Only {.arg year = 2010} is supported for now.")
    }
}



