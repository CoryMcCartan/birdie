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
           GEOID_zip = if_else(is.na(zip), str_c("cty", county), as.character(zip)))
# rm(voters)

p_r = prop.table(table(d$race))

geo_levels = c("county", "zip", "tract", "block")

r_probs = map(geo_levels, function(level) {
    predict_race_sgz(last_name, GEOID,
                     data=rename(d, GEOID=str_c("GEOID_", level)),
                     p_rgz=make_p_rgz(d, level), p_r=p_r)
}) %>%
    set_names(geo_levels)

# Party ID -----------

p_xr = prop.table(table(d$party, d$race))

fits = imap(r_probs, function(d_pr, level) {
    cat(level, "\n")

    lr = 0.4
    if (level == "block") level = "tract"
    if (level == "county") lr = 0.75
    if (level == "zip") lr = 0.5

    model_race(d_pr, party, !!rlang::sym(str_c("GEOID_", level)),
               data=d, config=list(lr=lr, tol_rhat=1.15, it_avgs=500))
})

write_rds(fits, here("data-out/nc_fits_party.rds"), compress="xz")



# Turnout -----------

d$n_voted = factor(d$n_voted)

p_xr = prop.table(table(d$n_voted, d$race))

fits = imap(r_probs, function(d_pr, level) {
    cat(level, "\n")

    lr = 0.65
    if (level == "block" || level == "tract") level = "zip"
    if (level == "county") lr = 0.8

    model_race(d_pr, n_voted, !!rlang::sym(str_c("GEOID_", level)),
               data=d, config=list(lr=lr, tol_rhat=1.2, draws=1000, it_avgs=500))
})

write_rds(fits, here("data-out/nc_fits_turnout.rds"), compress="xz")

