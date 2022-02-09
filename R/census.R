
# Helper function to make an R|S table
census_surname_table = function(S, S_name, p_r, regularize=TRUE, counts=FALSE) {
    if (length(p_r) != 5)
        cli_abort(c("Number of racial categories doesn't match the Census table.",
                    "i"="Categories should be White, Black, Hispanic, Asian, and other."))

    suppressMessages(require("wru"))
    x = data.frame(surname = unique(S)) %>%
        wru::merge_surnames(impute.missing=FALSE) %>%
        suppressWarnings() %>%
        na.omit

    x = rbind(x, list(surname="<generic>", surname.match="ALL OTHER NAMES",
                      p_whi=p_r[1], p_bla=p_r[2], p_his=p_r[3], p_asi=p_r[4], p_oth=p_r[5]))

    if (regularize | counts) {
        names_d = readRDS(system.file("data/names_2010_counts.rds", package="raceproxy"))
    }

    if (regularize) {
        x = left_join(x, names_d, by="surname.match") %>%
            mutate(count = dplyr::coalesce(count, 50))
        m = (0.2 + as.matrix(x[, 3:7]) * x$count) / (x$count + 1)
        x[, 3:7] = m
        x$count = NULL
    }
    if (counts) {
        x = left_join(x, names_d, by="surname.match") %>%
            mutate(count = dplyr::coalesce(count, 50))
        x[, 3:7] = as.matrix(x[, 3:7]) * (x$count + regularize)
        x$count = NULL
    }

    x = x[, -2]
    colnames(x) = c(S_name, "white", "black", "hisp", "asian", "other")
    as_tibble(x)
}

# Helper function to make an R|G table
census_zip_table = function(G, G_name, p_r, regularize=TRUE, count=FALSE) {
    if (!rlang::is_installed("zipWRUext2"))
        cli_abort(c("{.pkg zipWRUext2} must be installed to use ZIP code information automatically.",
                    ">"=' {.code devtools::install_github("https://github.com/jcuriel-unc/zipWRUext",
                    subdir="zipWRUext2")} to install.'))
    if (length(p_r) != 5)
        cli_abort(c("Number of racial categories doesn't match the Census data.",
                    "i"="Categories should be White, Black, Hispanic, Asian, and other."))

    x = zipWRUext2::zip_all_census2[, 2:8]
    x = rbind(x, list(zcta5="<none>", total_pop=1e4,
                      q_whi=1e4*p_r[1], q_bla=1e4*p_r[2],
                      q_his=1e4*p_r[3], q_asi=1e4*p_r[4], q_oth=1e4*p_r[5]))
    match_idx = match(levels(G), x$zcta5)
    match_idx[is.na(match_idx)] = match("<none>", x$zcta5)
    x = x[match_idx, ]
    alpha = if (regularize) c(3.1, 0.6, 0.8, 0.3, 0.2) else rep(0, 5)
    for (i in seq_along(p_r)) {
        denom = if (count) 1 else (x[, 2] + sum(alpha))
        x[, 2+i] = (x[, 2+i] + alpha[i]) / denom
    }
    x = x[, -2]
    colnames(x) = c(G_name, "white", "black", "hisp", "asian", "other")
    as_tibble(x)
}

