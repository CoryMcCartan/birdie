suppressMessages(library(here))
devtools::load_all(here("."))
library(wru)

if (file.exists(voterfile <- here("data/nc_voters.rds"))) {
    voters = readRDS(voterfile)
} else {
    voters = make_nc_statewide(voterfile)
}

load("data/censusData.rData")
censusData = censusData$NC
censusData$year = 2010
censusData = list(NC=censusData)

d_fit = slice_sample(voters, n=10e3) %>%
    mutate(gender = as.factor(coalesce(gender, "U"))) %>%
    filter(!is.na(age))
rm(voters); gc();
d_fit = d_fit %>%
    mutate(state = "NC",
           county = paste(county, "County")) %>%
    left_join(select(tigris::fips_codes, state, county, county_code), by=c("state", "county"))

# tables
p_xr = prop.table(table(d_fit$party, d_fit$race))
p_x = rowSums(p_xr)
p_r = colSums(p_xr)

r_probs_naive = d_fit %>%
    select(last=last_name, state, county=county_code) %>%
    predict_race_new(census.geo="county",
                     census.data=censusData) %>%
    as_tibble()

r_fit = select(r_probs_naive, starts_with("pred_"))
for (i in 1:5) r_fit[[i]][is.na(r_fit[[i]])] = p_r[i]
r_fit = as.matrix(r_fit) %>%
    apply(1, \(x) which.max(x) - 1L)

r_probs_me = d_fit %>%
    select(last=last_name, state, county=county_code) %>%
    mutate(across(everything(), as.character)) %>%
    predict_race_me("last",
                    census.geo="county",
                    #race.init = r_fit,
                    #census.surname=TRUE,
                    census.data=censusData,
                    control = list(iter = 500,
                                   burnin = 100,
                                   me.correct = "races")
                    ) %>%
    as_tibble()

