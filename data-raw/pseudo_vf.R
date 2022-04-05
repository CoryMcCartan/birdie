library(dplyr)
library(here)

if (file.exists(voterfile <- here("data-raw/nc_voters.rds"))) {
    voters = readRDS(voterfile)
} else {
    source(here("data-raw/nc.R"))
    voters = make_nc_statewide(voterfile)
}

r_probs = c(white=0.7, black=0.65, hisp=0.5, asian=0.5, aian=0.3, other=0.45)

pseudo_vf = voters %>%
    select(last_name, zip, race) %>%
    slice_sample(n=5000) %>%
    group_by(race) %>%
    mutate(zip = sample(zip),
           turnout = rbinom(n(), 1, r_probs[race[1]]),
           turnout = as.factor(if_else(turnout==1, "yes", "no"))) %>%
    ungroup() %>%
    mutate(across(where(is.factor), droplevels))

usethis::use_data(pseudo_vf, overwrite=TRUE, compress="xz")
