# Setup --------
FIPS_MIAMI = "086"
voters = read_rds(here("data-raw/l2_miami.rds")) |>
    mutate(tract = str_pad(tract, 6, pad="0"),
           party = fct_relevel(party, "dem", "ind", "rep", "lib"))

d_cens = census_race_geo_table("tract", state="FL", county=FIPS_MIAMI,
                               year=2019, survey="acs5") |>
    separate(GEOID, c("county", "tract"), sep=5)
d_county = group_by(d_cens, county) %>%
    summarize(across(white:other, sum), .groups="drop") %>%
    mutate(tract = NA)
d_geo = tigris::tracts(state="FL", county="086", cb=TRUE, year=2019) |>
    select(GEOID, geometry)

# make a R|GZ table for predict_race_sgz
make_p_rgx = function(voters) {
    voters_geo = distinct(voters, county, tract)
    group_by(d_cens, county, tract) %>%
        summarize(across(white:other, sum, na.rm=TRUE), .groups="drop") %>%
        bind_rows(d_county) %>%
        right_join(voters_geo, by=c("county", "tract")) %>%
        mutate(GEOID = ifelse(is.na(tract), as.character(county),
                              str_c(county, tract)),
               .before="county") %>%
        select(-county, -tract) %>%
        mutate(across(white:other, ~ (0.01 + .) /
                          (white + black + hisp + asian + other + aian + 0.06)))
}

# Subsample a table and format geography and outcome variables
set.seed(5118)

d = voters |>
    semi_join(bind_rows(d_cens, d_county), by=c("county", "tract")) |>
    mutate(GEOID = if_else(is.na(tract), as.character(county),
                           str_c(county, tract))) |>
    slice_sample(n=20e4) |>
    left_join(make_p_rgx(voters), by="GEOID")
print(as.character(head(d$last_name))) # ensure seed is working

# do BISG -----
p_r = with(voters, prop.table(table(race)))
p_xr = with(voters, prop.table(table(party, race), 2))

r_probs = bisg(~ nm(last_name) + GEOID, data=d, p_r=p_r,
               p_rgx=make_p_rgx(voters))


# do BIRDiE -----

fit = birdie(r_probs, party ~ white + black + hisp + (1 | GEOID), data=d)

wtd = est_weighted(r_probs, party ~ GEOID, data=d)
p_xr_thr = prop.table(table(d$party, predict(r_probs)), 2)

colSums(abs(coef(fit) - p_xr)) * 0.5
colSums(abs(coef(wtd) - p_xr)) * 0.5
colSums(abs(p_xr_thr - p_xr)) * 0.5

sum(abs((coef(fit) - p_xr) %*% diag(p_r))) * 0.5
sum(abs((coef(wtd) - p_xr) %*% diag(p_r))) * 0.5
sum(abs((p_xr_thr - p_xr) %*% diag(p_r))) * 0.5

