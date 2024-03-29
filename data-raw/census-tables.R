#############
# Prepare Census tables
library(easycensus)
library(tidyverse)
library(here)
devtools::load_all(".")


# National race over time -------
years = 2005:2021
get_natl_race <- function(year) {
    if (year == 2010) {
        x = census_race_geo_table("us", year=2010, survey="dec", counts=FALSE)
    } else if (year == 2020) {
        return(NULL)
    } else {
        x = census_race_geo_table("us", year=year, survey="acs1", counts=FALSE)
    }

    out = as.numeric(x[, -1:-2])
    names(out) = c("white", "black", "hisp", "asian", "aian", "other")
    out
}

l_race_year = map(years, get_natl_race)
names(l_race_year) = years

l_race_year[["2020"]] = l_race_year[["2019"]]*0.5 + l_race_year[["2021"]]*0.5
# run sysdata.R to save!



# State by race ------
d = census_race_geo_table("state", year=2010, survey="dec", counts=TRUE, GEOIDs=TRUE)
write_rds(d, here("inst/extdata/state_race_2010.rds"), compress="xz")


# ZIP code by race ------
d = census_race_geo_table("zcta", year=2010, survey="dec", counts=TRUE, GEOIDs=TRUE)
write_rds(d, here("inst/extdata/zip_race_2010.rds"), compress="xz")


# ZIP code by race (2020) ------
d_raw = censable::build_dec("block", "NC", geometry=TRUE, year=2020)
zip_geom = tigris::zctas(state="NC", year=2010)
idx_zip = geomander::geo_match(d_raw, zip_geom, method="area")
d = d_raw |>
    sf::st_drop_geometry() |>
    mutate(zip = zip_geom$ZCTA5CE10[idx_zip],
           pop_asian = pop_asian + pop_nhpi,
           pop_other = pop_other + pop_two,
           vap_asian = vap_asian + vap_nhpi,
           vap_other = vap_other + vap_two) |>
    select(-pop_nhpi, -pop_two, -vap_nhpi, -vap_two) |>
    group_by(zip) |>
    summarize(across(pop:vap_other, function(x) as.integer(sum(x))))
write_rds(d, here("inst/extdata/zip_race_2020_nc.rds"), compress="xz")


# Surnames by race ------
# overall race probabilities
p_r = summarize(d, across(pop:pop_other, sum)) |>
    mutate(across(everything(), ~ . / pop)) |>
    select(-pop) |>
    rename_with(~ paste0("pr", str_sub(., 4))) |>
    as.matrix() |>
    `[`(1, )

if (!file.exists(path <- here("data-raw/names/Names_2010Census.csv"))) {
    url <- "https://www2.census.gov/topics/genealogy/2010surnames/names.zip"
    zipfile <- here("data-raw/names.zip")
    download.file(url, zipfile)
    dir.create(dirname(path), showWarnings=FALSE)
    unzip(zipfile, files="Names_2010Census.csv", exdir=dirname(path))
    file.remove(zipfile)
}
d_raw = read_csv(path, col_types="c-i--dddddd", na=c("", "NA", "(S)"))


d = d_raw |>
    select(last_name=name, count,
           pr_white=pctwhite,
           pr_black=pctblack,
           pr_hisp=pcthispanic,
           pr_asian=pctapi,
           pr_aian=pctaian,
           pr_other=pct2prace) |>
    mutate(across(pr_white:pr_other, ~ ./100))
m = select(d, starts_with("pr_")) |>
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
d = select(d, last_name, count) |>
    bind_cols(m)
d$last_name[nrow(d)] = "<generic>"

# convert to ints
d = d |>
    mutate(across(starts_with("pr_"), ~ as.integer(. * count))) |>
    select(-count)

write_rds(d, here("inst/extdata/names_2010_counts.rds"), compress="xz")
