
# Helper function to make an R|S table
census_surname_table = function(S, S_name, p_r, regularize=TRUE, counts=FALSE) {
    if (length(p_r) != 6)
        cli_abort(c("Number of racial categories doesn't match the Census table.",
                    "i"="Categories should be White, Black, Hispanic, Asian,
                    American Indian/Alaska Native, and other."))

    d_cens = readRDS(system.file("extdata", "names_2010_counts.rds",
                                 package="raceproxy", mustWork=TRUE))

    out = data.frame(last_name = unique(proc_names(S))) %>%
        left_join(d_cens, by="last_name")
    missing_idx = which(is.na(out$pr_white))
    double_idx = missing_idx[is_double_name(out$last_name[missing_idx])]
    # track which indices will become <generic>
    bad = setdiff(missing_idx, double_idx)

    # double-barelled surnames
    match_1st_idx = match(stringr::word(out$last_name[double_idx], 1, 1), d_cens$last_name)
    match_2nd_idx = match(stringr::word(out$last_name[double_idx], 2, 2), d_cens$last_name)

    # neither matches
    na_ct = is.na(match_1st_idx) + is.na(match_2nd_idx)
    bad = c(bad, double_idx[which(na_ct == 2)])

    # one or the other matches: use it but with half confidence
    only_1st = which(!is.na(match_1st_idx) & is.na(match_2nd_idx))
    only_2nd = which(is.na(match_1st_idx) & !is.na(match_2nd_idx))
    out[double_idx[only_1st], -1] = d_cens[match_1st_idx[only_1st], -1] / 2
    out[double_idx[only_2nd], -1] = d_cens[match_2nd_idx[only_2nd], -1] / 2

    # both match: average guesses
    both = which(na_ct == 0)
    pr_1st = as.matrix(d_cens[match_1st_idx[both], -1])
    pr_2nd = as.matrix(d_cens[match_2nd_idx[both], -1])
    sum1 = rowSums(pr_1st)
    sum2 = rowSums(pr_2nd)
    new_totals = pmin(sum1, sum2) / 2 # err on the side of less information
    out[double_idx[both], -1] =  0.5*new_totals*(pr_1st/sum1 + pr_2nd/sum2)

    # replace non-matches w/ <generic>
    out = rbind(out[-bad, ], tail(d_cens, 1))

    if (regularize) {
        out[, -1] = 0.2 + out[, -1]
    }
    if (!counts) {
        out[, -1] = out[, -1] / (rowSums(out[, -1]) + 0.2*regularize)
    }

    colnames(out) = c(S_name, "white", "black", "hisp", "asian", "aian", "other")
    as_tibble(out)
}

# Helper function to make an R|G table
census_zip_table = function(G, G_name, p_r, regularize=TRUE, counts=FALSE) {
    if (length(p_r) != 6)
        cli_abort(c("Number of racial categories doesn't match the Census table.",
                    "i"="Categories should be White, Black, Hispanic, Asian,
                    American Indian/Alaska Native, and other."))

    d_cens = readRDS(system.file("extdata", "zip_race_2010.rds",
                            package="raceproxy", mustWork=TRUE))
    d_cens = rbind(d_cens, list(zcta5="<none>", pop=1e4,
                      pop_white=1e4*p_r[1], pop_black=1e4*p_r[2],
                      pop_hisp=1e4*p_r[3], pop_asian=1e4*p_r[4],
                      pop_aian=1e4*p_r[5], pop_other=1e4*p_r[6]))
    match_idx = match(levels(G), d_cens$zcta5)
    match_idx[is.na(match_idx)] = match("<none>", d_cens$zcta5)
    d_cens = d_cens[match_idx, ]
    alpha = if (regularize) c(3.1, 0.60, 0.87, 0.24, 0.036, 0.11) else rep(0, 6)
    for (i in seq_along(p_r)) {
        denom = if (counts) 1 else (d_cens[, 2] + sum(alpha))
        d_cens[, 2+i] = (d_cens[, 2+i] + alpha[i]) / denom
    }
    d_cens = d_cens[, -2]
    colnames(d_cens) = c(G_name, "white", "black", "hisp", "asian", "aian", "other")
    as_tibble(d_cens)
}

