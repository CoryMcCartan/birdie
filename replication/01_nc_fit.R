voters = read_rds(here("data-raw/nc_voters.rds"))
d_cens = read_rds(here("data-raw/nc_block_vap_race.rds")) %>%
    drop_na()

# make a R|GZ table for predict_race_sgz
make_p_rgz = function(voters, level=c("block", "tract", "county", "zip"), counts=FALSE) {
    d_county = group_by(d_cens, county) %>%
        summarize(across(vap:vap_other, sum), .groups="drop") %>%
        mutate(tract = NA, block = NA)

    if (level == "block") {
        voters_geo = distinct(voters, county, tract, block)
        d = bind_rows(d_cens, d_county) %>%
            right_join(voters_geo, by=c("county", "tract", "block")) %>%
            mutate(GEOID = ifelse(is.na(block), county, str_c(county, tract, block)),
                   .before="county") %>%
            select(-county, -tract, -block)
    } else if (level == "tract") {
        voters_geo = distinct(voters, county, tract)
        d = group_by(d_cens, county, tract) %>%
            summarize(across(vap:vap_other, sum), .groups="drop") %>%
            bind_rows(d_county) %>%
            right_join(voters_geo, by=c("county", "tract")) %>%
            mutate(GEOID = ifelse(is.na(tract), county, str_c(county, tract)),
                   .before="county") %>%
            select(-county, -tract, -block)
    } else if (level == "zip") {
        voters_geo = voters %>%
            distinct(county, zip) %>%
            mutate(GEOID = if_else(is.na(zip), str_c("cty", county), as.character(zip))) %>%
            distinct(GEOID)
        d_county = mutate(d_county, GEOID = str_c("cty", county)) %>%
            select(-county, -tract, -block)
        d = census_zip_table(NULL, "GEOID", 1:6, counts=TRUE) %>%
            rename_with(~ str_c("vap_", .), white:other) %>%
            mutate(vap = rowSums(across(vap_white:vap_other))) %>%
            bind_rows(d_county) %>%
            right_join(voters_geo, by="GEOID")
    } else {
        d = rename(d_county, GEOID = county) %>%
            select(-tract, -block)
    }

    if (isFALSE(counts)) {
        d = mutate(d, across(vap_white:vap_other, ~ (0.01 + .) / (vap + 0.06)))
    }

    d %>%
        select(-vap) %>%
        rename_with(~ str_sub(., 5), vap_white:vap_other) %>%
        mutate(GEOID = as.character(GEOID))
}

set.seed(5118)
d = slice_sample(voters, n=500e3) %>%
    mutate(GEOID_county = as.character(county),
           GEOID_tract = if_else(is.na(tract), GEOID_county, str_c(county, tract)),
           GEOID_block = if_else(is.na(block), GEOID_county, str_c(county, tract, block)),
           GEOID_zip = if_else(is.na(zip), str_c("cty", county), as.character(zip)),
           n_voted = factor(n_voted),
           party=coalesce(party, "ind"))
rm(voters)

p_r = prop.table(table(d$race))

geo_levels = c("county", "zip", "tract", "block")

r_probs = map(geo_levels, function(level) {
    predict_race_sgz(last_name, GEOID,
                     data=rename(d, GEOID=str_c("GEOID_", level)),
                     p_rgz=make_p_rgz(d, level), p_r=p_r)
}) %>%
    set_names(geo_levels)

# BISG quality
log_score_baseline = mean(log(p_r[as.integer(d$race)]))
log_scores = map_dbl(r_probs, function(x) {
    pr_act = as.matrix(x)[cbind(1:nrow(x), as.integer(d$race))]
    pr_act[pr_act == 0] = 1e-6
    mean(log(pr_act))
})
acc_thresh = map_dbl(r_probs, ~ mean(max.col(.) == as.integer(d$race)))

list(base_score = log_score_baseline,
     score = log_scores,
     acc = acc_thresh) |>
    write_rds(here("paper/data/nc_bisg.rds"))


# Party ID -----------

if (!file.exists(path <- here("data-out/nc_fits_party.rds"))) {
    fits_party = list()
    # fits_party = imap(r_probs, function(d_pr, level) {
    for (level in names(r_probs)) {
        d_pr = r_probs[[level]]
        cat(level, "\n")

        lr = 0.4
        if (level == "block") level = "tract"
        if (level == "county") lr = 0.75
        if (level == "zip") lr = 0.5

        # fits_party[[level]] = model_race(d_pr, party, !!rlang::sym(str_c("GEOID_", level)),
        fits_party$block = model_race(d_pr, party, !!rlang::sym(str_c("GEOID_", level)),
                                      data=d, config=list(lr=lr, tol_rhat=1.15, it_avgs=500))
    # })
        gc()
    }

    write_rds(fits_party, path, compress="xz")
} else {
    fits_party <- read_rds(path)
}


# Turnout -----------

pre12 = which(d$reg_date <= as.Date("2012-11-01"))
if (!file.exists(path <- here("data-out/nc_fits_turnout.rds"))) {
    fits_turnout = imap(r_probs, function(d_pr, level) {
        cat(level, "\n")

        lr = 0.65
        if (level == "block" || level == "tract") level = "zip"
        if (level == "county") lr = 0.8

        model_race(d_pr[pre12, ], n_voted, !!rlang::sym(str_c("GEOID_", level)),
                   data=d[pre12, ], config=list(lr=lr, tol_rhat=1.2, draws=1000, it_avgs=500))
    })

    write_rds(fits_turnout, path, compress="xz")
} else {
    fits_turnout <- read_rds(path)
}

