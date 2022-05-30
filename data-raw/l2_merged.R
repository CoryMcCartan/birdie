library(here)
library(dplyr)
library(stringr)

source(here("data-raw/nc.R"))

ncvf = make_nc_statewide(NULL)
l2vf = readr::read_rds(here("data-raw/VF/l2vf_nc.rds"))

ncvf = distinct(ncvf, last_name, first_name, county, race, .keep_all=TRUE) %>%
    select(-regnum)
l2vf = distinct(l2vf, last_name, first_name, county, race, .keep_all=TRUE) %>%
    rename(zip_l2=zip,
           gender_l2=gender)

l2vf = l2vf %>%
    mutate(first_name = str_to_upper(first_name),
           last_name = str_to_upper(last_name))

d = left_join(ncvf, l2vf, by=c("county", "race", "last_name", "first_name"))
d$gender[is.na(d$gender)] = "U"

d %>%
    mutate(across(where(is.character), factor)) %>%
    relocate(tract:block, .after=county) %>%
    readr::write_rds(here("data-raw/nc_voters.rds"), compress="xz")

rm(ncvf, l2vf)

d_cens = censable::build_dec("block", "NC", year=2010, geometry=FALSE)
d_cens %>%
    mutate(county = factor(str_sub(GEOID, 1, 5), levels=levels(d$county)),
           tract = factor(str_sub(GEOID, 6, 11), levels=levels(d$tract)),
           block = factor(str_sub(GEOID, 12, 15), levels=levels(d$block)),
           vap_asian = vap_asian + vap_nhpi,
           vap_other = vap_other + vap_two,
           across(vap:vap_other, as.integer)) %>%
    select(county:block, vap, vap_white, vap_black, vap_hisp, vap_asian, vap_aian, vap_other) %>%
    readr::write_rds(here("data-raw/nc_block_vap_race.rds"), compress="xz")
