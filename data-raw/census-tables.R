#############
# Prepare Census tables
library(easycensus)
library(tidyverse)
library(here)


# ZIP code by race ------
# find_dec_table("race")
d_raw = get_dec_table("zcta", "P005")

d = d_raw %>%
    transmute(zcta5 = GEOID,
              value = value,
              race = case_when(
                  race == "total" ~ "total",
                  hispanic_or_latino_origin == "hispanic or latino" ~ "hisp",
                  TRUE ~ as.character(tidy_race(race))
                  ),
              race = fct_collapse(race,
                                  asian=c("asian", "nhpi"),
                                  other=c("other", "two"))) %>%
    group_by(zcta5, race) %>%
    summarize(value = sum(value),
              .groups="drop") %>%
    pivot_wider(names_from=race) %>%
    select(zcta5, pop=total, pop_white=white, pop_black=black, pop_hisp=hisp,
           pop_asian=asian, pop_aian=aian, pop_other=other) %>%
    mutate(across(-zcta5, as.integer))

write_rds(d, here("inst/extdata/zip_race_2010.rds"), compress="xz")


# ZIP code by race (2020) ------
d_raw = censable::build_dec("block", "NC", geometry=TRUE, year=2020)
zip_geom = tigris::zctas(state="NC", year=2010)
idx_zip = geomander::geo_match(d_raw, zip_geom, method="area")
d = d_raw %>%
    sf::st_drop_geometry() %>%
    mutate(zip = zip_geom$ZCTA5CE10[idx_zip],
           pop_asian = pop_asian + pop_nhpi,
           pop_other = pop_other + pop_two,
           vap_asian = vap_asian + vap_nhpi,
           vap_other = vap_other + vap_two) %>%
    select(-pop_nhpi, -pop_two, -vap_nhpi, -vap_two) %>%
    group_by(zip) %>%
    summarize(across(pop:vap_other, function(x) as.integer(sum(x))))
write_rds(d, here("inst/extdata/zip_race_2020_nc.rds"), compress="xz")


# Surnames by race ------
# overall race probabilities
p_r = summarize(d, across(pop:pop_other, sum)) %>%
    mutate(across(everything(), ~ . / pop)) %>%
    select(-pop) %>%
    rename_with(~ paste0("pr", str_sub(., 4))) %>%
    as.matrix() %>%
    `[`(1, )

d_raw = read_csv(here("data-raw/names/Names_2010Census.csv"),
                 col_types="c-i--dddddd", na=c("", "NA", "(S)"))

d = d_raw %>%
    select(last_name=name, count,
           pr_white=pctwhite,
           pr_black=pctblack,
           pr_hisp=pcthispanic,
           pr_asian=pctapi,
           pr_aian=pctaian,
           pr_other=pct2prace) %>%
    mutate(across(pr_white:pr_other, ~ ./100))
m = select(d, starts_with("pr_")) %>%
    as.matrix()

# allocate missing row total to NAs in proportion to population totals
rowtot = rowSums(m, na.rm=T)
rowwgt = as.numeric(is.na(m) %*% p_r)
for (i in seq_len(ncol(m))) {
    m[, i] = coalesce(m[, i], (p_r[i] / rowwgt) * (1 - rowtot))
}
# ensure sum to 1
rowtot = rowSums(m, na.rm=T)
m = m / rowtot
# format as tibble and clean up
colnames(m) = names(p_r)
d = select(d, last_name, count) %>%
    bind_cols(m)
d$last_name[nrow(d)] = "<generic>"

# convert to ints
d = d %>%
    mutate(across(starts_with("pr_"), ~ as.integer(. * count))) %>%
    select(-count)

write_rds(d, here("inst/extdata/names_2010_counts.rds"), compress="xz")


# national tract table
d_raw = purrr::map_dfr(state.abb, ~ censable::build_dec("tract", ., geometry=FALSE, year=2010))

d_raw %>%
    mutate(pop_asian = pop_asian + pop_nhpi,
           pop_other = pop_other + pop_two,
           vap_asian = vap_asian + vap_nhpi,
           vap_other = vap_other + vap_two) %>%
    select(tract=GEOID, white=vap_white, black=vap_black, hisp=vap_hisp,
           asian=vap_asian, aian=vap_aian, other=vap_other) %>%
    mutate(across(white:other, as.integer)) %>%
    readr::write_csv("~/Desktop/tract_race_2010.csv")

