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
#'   can be identified similarly: `zip()` for ZIP codes, and `state()` for
#'   states. If no other predictor variables are provided, then `bisg()` will
#'   automatically be able to build a table of census data to use in inference.
#'   If other predictor variables are included, or if other geographic
#'   identifiers are used, then the user must specify the `p_rgx` argument below.
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
#'   The [census_race_geo_table()] function can be helpful to prepare tables,
#'   as can be the `build_dec()` and `build_acs()` functions in the `censable`
#'   package.
#' @param p_rs The distribution of race given last name. As with `p_rgx`, should
#'   be provided as a data frame, with a column of names and additional columns
#'   for each racial group. Users should not have to specify this argument in
#'   most cases, as the table will be built from published Census surname tables
#'   automatically. Counts are required for `bisg_me()`.
#' @param save_rgx If `TRUE`, save the `p_rgx` table (matched to each
#'   individual) as the `"p_rgx"` and `"gx"` attributes of the output.
#'   Necessary for some sensitivity analyses.
#'
#' @return An object of class `bisg`, which is just a data frame with some
#'   additional attributes. The data frame has rows matching the input data and
#'   columns for the race probabilities.
#'
#' @examples
#' data(pseudo_vf)
#' bisg(~ nm(last_name), data=pseudo_vf)
#'
#' r_probs = bisg(~ nm(last_name) + zip(zip), data=pseudo_vf)
#' summary(r_probs)
#' head(predict(r_probs))
#'
#' @references
#' Elliott, M. N., Fremont, A., Morrison, P. A., Pantoja, P., and Lurie, N.
#' (2008). A new method for estimating race/ethnicity and associated disparities
#' where administrative records lack self-reported race/ethnicity. *Health
#' Services Research*, 43(5p1):1722–1736.
#'
#' Fiscella, K. and Fremont, A. M. (2006). Use of geocoding and surname analysis
#' to estimate race and ethnicity. *Health Services Research*,
#' 41(4p1):1482–1500.
#'
#' Imai, K., Olivella, S., & Rosenman, E. T. (2022). Addressing census data
#' problems in race imputation via fully Bayesian Improved Surname Geocoding and
#' name supplements. *Science Advances*, 8(49), eadc9824.
#'
#' @describeIn bisg The standard BISG model.
#' @concept bisg
#' @export
bisg <- function(formula, data=NULL, p_r=p_r_natl(), p_rgx=NULL, p_rs=NULL,
                 save_rgx=TRUE) {
    vars = parse_bisg_form(formula, data)

    l_name = make_name_tbl_vec(vars, p_r, p_rs, FALSE)
    l_gx = make_gx_tbl_vec(vars, p_r, p_rgx)

    m_bisg = est_bisg(l_name$S, l_gx$GX, l_name$p_sr, l_gx$p_gxr, l_gx$p_r)

    out = as_tibble(m_bisg)
    class(out) = c("bisg", class(out))
    attr(out, "S_name") = vars$S_name
    attr(out, "GX_names") = colnames(vars$GX)
    attr(out, "p_r") = l_gx$p_r
    if (isTRUE(save_rgx)) {
        attr(out, "p_rgx") = l_gx$p_rgx
        attr(out, "gx") = l_gx$GX
    }
    attr(out, "method") = "std"

    out
}

#' @param iter How many sampling iterations in the Gibbs sampler
#' @param warmup How many burn-in iterations in the Gibbs sampler
#' @param cores How many parallel cores to use in computation. Around 4 seems to
#'   be optimal, even if more are available.
#'
#' @examples
#' data(pseudo_vf)
#' bisg_me(~ nm(last_name) + zip(zip), data=pseudo_vf)

#' @describeIn bisg The measurement error BISG model.
#' @concept bisg
#' @export
bisg_me <- function(formula, data=NULL, p_r=p_r_natl(), p_rgx=NULL, p_rs=NULL,
                    iter=1000, warmup=100, cores=1L) {
    vars = parse_bisg_form(formula, data)

    l_name = make_name_tbl_vec(vars, p_r, p_rs, TRUE)
    l_gx = make_gx_tbl_vec(vars, p_r, p_rgx)

    # set up sampler vectors
    n_gx = nrow(l_gx$p_rgx)
    n_s = nrow(l_name$p_sr)
    alpha_gzr = matrix(rep(l_gx$p_r, n_gx), nrow=n_gx, ncol=6, byrow=T)
    beta_sr = matrix(rep(l_gx$p_r, n_s), nrow=n_s, ncol=6, byrow=T)

    # computation params
    iter = as.integer(max(iter, 1))
    warmup = as.integer(max(warmup, 1))
    cores = if (cores == 1) 0L else as.integer(cores)

    # extra process_done to try to handle cases with interruptions
    cli::cli_process_done()
    m_bisg = gibbs_me(iter+warmup, warmup, l_name$S, l_gx$GX,
                      l_name$p_sr, l_gx$p_rgx,
                      alpha_gzr, beta_sr, cores=cores, verbosity=3L)
    cli::cli_process_done()
    colnames(m_bisg) = paste0("pr_", names(p_r))

    out = as_tibble(m_bisg)
    class(out) = c("bisg_me", "bisg", class(out))
    attr(out, "S_name") = vars$S_name
    attr(out, "GX_names") = colnames(vars$GX)
    attr(out, "p_r") = l_gx$p_r
    attr(out, "method") = "me"

    out
}



# Parse formula and return a list with name vector, name vector name,
# type of geography used (or "none"), and data frame of other covariates
parse_bisg_form <- function(formula, data=NULL) {
    formula = update.formula(formula, NULL ~ .)
    f_terms = terms(formula, specials=c("nm", "zip", "state"), data=data)
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
    if (!check_discrete(S)) {
        cli_abort("The names vector provided in {.fn nm} must be
                   a character or factor with no missing values.",
                  call=parent.frame())
    }

    # check ZIP, state, etc.
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
        # process names (by unique level)
        idx_uniq = vctrs::vec_unique_loc(S)
        idx_dup = to_unique_ids(S)
        S = proc_name(S[idx_uniq], to_latin=TRUE)[idx_dup]

        if (is.character(p_r)) p_r = p_r_natl()
        if (length(p_r) != 6) {
            cli_abort(c("Number of racial categories doesn't match the Census table.",
                        "i"="Categories should be `white`, `black`, `hisp`, `asian`,
                        `aian` (American Indian/Alaska Native), and `other`."),
                      call=parent.frame())
        }

        if (isTRUE(for_me)) {
            p_rs = census_surname_table(S, vars$S_name, counts=TRUE, flip=FALSE)
        } else {
            p_rs = census_surname_table(S, vars$S_name, counts=FALSE, flip=TRUE)
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
        if (p_r != "estimate" && length(p_r) != ncol(p_rs) + 1) {
            cli_abort("Number of racial categories in {.arg p_rs}
                      and {.arg p_r} must match.", call=parent.frame())
        }

        S = factor(S, levels=p_rs[[name_col]])

        p_s = prop.table(table(S))
        p_rs = as.matrix(p_rs[, -name_col])
        p_sr = p_rs / rowSums(p_rs)
        for (i in seq_len(ncol(p_sr))) {
            p_sr[, i] = p_sr[, i] * p_s
            p_sr[, i] = p_sr[, i] / sum(p_sr[, i])
        }
    }

    # reorder columns to match `p_r`
    if (p_r[1] != "estimate") {
        if (is.null(names(p_r))) cli_abort("{.arg p_r} must have names.", call=parent.frame())
        if (is.null(colnames(p_sr))) cli_abort("{.arg p_rs} must have column names.", call=parent.frame())
        idx_names = match(names(p_r), colnames(p_sr))
        if (any(is.na(idx_names))) {
            cli_abort(c("Names of {.arg p_r} and column names of {.arg p_rs} must match.",
                        "i"="Names for {.arg p_r}: {.val {names(p_r)}}",
                        "i"="Names for {.arg p_rs}: {.val {names(p_rs)}}",
                        ">"="If you provided {.arg p_r} but not {.arg p_rs}, make sure
                    {.arg p_r} has {.fn names} matching `white`, `black`,
                    `hisp`, `asian`, `aian`, and `other`."),
                    call=parent.frame())
        }
        p_sr = p_sr[, idx_names]
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

    est_p_r = is.character(p_r) && (p_r == "est" || p_r == "estimate")

    p_r_tmp = p_r
    if (est_p_r) p_r_tmp = p_r_natl()
    GX_names = colnames(vars$GX)
    if (vars$geo_type != "none") {
        if (length(p_r_tmp) != 6) {
            cli_abort(c("Number of racial categories doesn't match the Census table.",
                        "i"="Categories should be White, Black, Hispanic, Asian,
                    American Indian/Alaska Native, and other."), call=parent.frame())
        }

        if (vars$geo_type == "zip") {
            p_rgx = census_premade_table(GX_names, "zip_race_2010.rds", counts=TRUE)

            # match ZIPs to ZCTAs
            vars$GX[[1]] = proc_zip(vars$GX[[1]])
        } else if (vars$geo_type == "state") {
            p_rgx = census_premade_table(GX_names, "state_race_2010.rds", counts=TRUE)

            # match states
            vars$GX[[1]] = proc_state(vars$GX[[1]])
        } else {
            cli_abort(c("{.fn {vars$geo_type}} is not yet supported.",
                        ">"="Please file an issue at
                        {.url https://github.com/CoryMcCartan/birdie/issues/new}"),
                        call=parent.frame())
        }

    } else { # "none" -> use `p_rgx`
        if (!is.data.frame(p_rgx)) {
            cli_abort("{.arg p_rgx} must be a data frame.", call=parent.frame())
        }
        if (!all(colnames(vars$GX) %in% colnames(p_rgx))) {
            cli_abort("All predictor columns must be in {.arg p_rgx}.", call=parent.frame())
        }
        if (!is.character(p_r) && ncol(p_rgx) != length(p_r) + ncol(vars$GX)) {
            cli_abort("Number of racial categories in {.arg p_rgx}
                      and {.arg p_r} must match.", call=parent.frame())
        }
        if (anyDuplicated(p_rgx) > 0) {
            cli_abort("{.arg p_rgx} must have unique rows.", call=parent.frame())
        }

        # if there is just one ID column we can easily fill in "<none>"
        if (length(GX_names) > 1) {
            d_miss <- distinct(anti_join(vars$GX, p_rgx, by=names(vars$GX)))
            if (nrow(d_miss) > 0) {
                str_miss = capture.output(head(d_miss, 10))
                if (nrow(d_miss) > 10) str_miss = c(str_miss, "  ...")
                msg = c("Some predictor combinations are missing from {.arg p_rgx}:", str_miss)
                mgs = str_replace_all(msg, " ", "\ua0")
                rlang::abort(msg, use_cli_format=TRUE, call=parent.frame())
            }
        }
    }

    # tidy up
    if (length(GX_names) == 1) {
        p_rgx[[GX_names]] = as.character(p_rgx[[GX_names]])
        vars$GX[[1]] = as.character(vars$GX[[1]])

        if (!"<none>" %in% p_rgx[[GX_names]]) {
            new_row = as.data.frame(as.list(1e5 * p_r_tmp))
            new_row[[GX_names]] = "<none>"
            p_rgx = rbind(p_rgx, new_row)
        }

        rowtot = rep_len(0, nrow(p_rgx))
        for (nm in names(p_r_tmp)) {
            p_rgx[[nm]] = coalesce(p_rgx[[nm]], p_r_tmp[nm])
            rowtot = rowtot + p_rgx[[nm]]
        }
        p_rgx = p_rgx[rowtot > 0, ] # drop rows with no people

        match_idx = match(vars$GX[[1]], p_rgx[[GX_names]])
        idx_miss = which(is.na(match_idx))
        if (length(idx_miss) > 0) {
            vars$GX[[1]][idx_miss] = "<none>"
        }
    }

    # get columns in order
    idx_names = setdiff(seq_len(ncol(p_rgx)), match(GX_names, colnames(p_rgx))) # default to non-GX-names cols
    if (!est_p_r) {
        idx_names = match(names(p_r), colnames(p_rgx))
        if (any(is.na(idx_names))) {
            cli_abort(c("Names of {.arg p_r} and column names of {.arg p_rgx} must match.",
                        "i"="Names for {.arg p_r}: {.val {names(p_r)}}",
                        "i"="Names for {.arg p_rgx}: {.val {names(p_rgx)}}",
                        ">"="If you provided {.arg p_r} but not {.arg p_rgx}, make sure
                        {.arg p_r} has {.fn names} matching `white`, `black`,
                        `hisp`, `asian`, `aian`, and `other`."),
                      call=parent.frame())
        }
    }

    # subset to needed rows
    d_match = left_join(vars$GX, p_rgx, by=GX_names)
    GX_vec = as.factor(vctrs::vec_duplicate_id(d_match))
    p_rgx = as.matrix(d_match[vctrs::vec_unique_loc(GX_vec), idx_names])

    # flip which margin sums to 1
    p_gx = prop.table(table(GX_vec))
    p_gxr = p_rgx / rowSums(p_rgx)
    for (i in seq_len(ncol(p_gxr))) {
        p_gxr[, i] = p_gxr[, i] * p_gx
        p_gxr[, i] = p_gxr[, i] / sum(p_gxr[, i])
    }

    if (est_p_r) {
        p_r = colSums(p_rgx)
        p_r = p_r / sum(p_r)
    }

    list(GX = GX_vec,
         p_r = p_r,
         p_rgx = p_rgx,
         p_gxr = p_gxr)
}

# Call Bayes' rule C++
est_bisg = function(S, GX, p_sr, p_gxr, p_r, geo=TRUE) {
    if (!geo) p_gxr = p_gxr*0 + 1

    m_bisg = calc_bayes_bisg(S, GX, p_sr, p_gxr, p_r);
    colnames(m_bisg) = paste0("pr_", names(p_r))
    m_bisg
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
p_r_natl <- function(year=2021, vap=FALSE) {
    if (isTRUE(vap)) cli_abort("Only {.arg vap = FALSE} is supported for now.")

    year = as.integer(year)
    if (year < 2005L || year > 2021L) {
        cli_abort("{.arg year} must be between 2005 and 2021 (inclusive).")
    }

    l_race_year[[as.character(year)]]
}
