library(birdie)
library(readr)

# STEP 1: do BISG -------------------------------------------------------------

## Option 1: produce BISG estimates in R -----

# a data frame with columns containing
# - outcome variable (HMID indicator). Assume this is called `hmid` in the rest of the template
# - last name. Assume we call this `last_name`
# - BISG geo unit, e.g. tract. Assume we call this `GEOID`
# - model geo unit. probably shouldn't be as granular as block group -- try counties?
#   truncate tract FIPS code with `county = str_sub(GEOID, 1, 5)
#   Assume this is named `county`
#   No missing values are allowed here. So those need to be filled in with a
#   state code or (worst case) an "Unknown" level for the model to run properly.
itm_data = ...

# a data frame describing the racial breakdown of each geo unit
# - first column should be named `GEOID` (i.e., it should match the name of the
#   geo column in `itm_data`)
# - other columns should be named `white`, `black`, `hisp`, `asian`, `aian`, and `other`
# - can use `census_race_geo_table()` to make this: e.g. this will get national tract-level data from 2021 ACS
p_rg_table = do.call(rbind, lapply(state.abb, function(state) {
    census_race_geo_table("tract", state=state, year=2021, survey="acs5")
}))

# vector of prior probability estimates for the overall tax-filing population
# these are taken from the March 2021 CPS ASEC
p_r = c(white=0.633, black=0.118, hisp=0.165, asian=0.063, aian=0.00745, other=0.014)

# substitute `bisg_me()` for measurement-error-improved probabilities
r_probs = bisg(~ nm(last_name) + GEOID, data=itm_data,
               p_r=p_r, p_rgx=p_rg_table)

## Option 2: load in BISG probabilities from CSV ------
# columns should be `pr_white`, `pr_black`, `pr_hisp`, `pr_asian`, `pr_aian`, and `pr_other`
r_probs = read_csv(...)


# STEP 2: Fit the model  -------------------------------------------------------

# start with this
fit = birdie(r_probs, hmid ~ county, data=itm_data)

# if too slow, you can pass in `ctrl = birdie.ctrl(abstol=1e-5)`, which will increase numerical tolerances slightly

# try this too (preferred, if it's not too slow)
fit = birdie(r_probs, hmid ~ (1 | county), data=itm_data)

# if you have any county-level predictors (% of county that are taxpayers,
# socioeconomic variables, anything else), you can include them & it will help
# produce better estimates:
fit = birdie(r_probs, hmid ~ pct_taxpaying + (1 | county), data=itm_data)

# get some idea of the standard errors
birdie(r_probs, hmid ~ 1, data=itm_data, se_boot=150)$se

# STEP 3: Summarize estimates  -------------------------------------------------
library(dplyr)

coef(fit) # show overall estimates of outcome by race
plot(fit) # visualize the estimates

# see all pairwise racial disparities
disparities(fit) %>%
    filter(disparity > 0) %>%
    arrange(desc(disparity))

# produce state-level estimates
r_probs_y = fitted(fit) # updated BISG probabilities (conditional on outcome as well)
est_weighted(r_probs_y, party ~ state, data=itm_data) %>%
    tidy(subgroup = TRUE)

# can also get ZIP or tract-level estimates with `tidy(fit, subgroup=TRUE)`,
# and then aggregate those to the state-level
