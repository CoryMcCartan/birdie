
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

