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
    mutate(id = row_number()) %>%
    as.data.frame() %>%
    predict_race_new(census.geo="county",
                     census.data=censusData) %>%
    as_tibble() %>%
    arrange(id)

r_fit = select(r_probs_naive, starts_with("pred_"))
for (i in 1:5) r_fit[[i]][is.na(r_fit[[i]])] = p_r[i]
r_fit = as.matrix(r_fit) %>%
    apply(1, \(x) which.max(x) - 1L)

r_probs_me = d_fit %>%
    select(last=last_name, state, county=county_code) %>%
    mutate(across(everything(), as.character),
           id = row_number()) %>%
    as.data.frame() %>%
    predict_race_me("last",
                    census.geo="county",
                    race.init = r_fit,
                    census.data=censusData,
                    control = list(iter = 500,
                                   burnin = 100,
                                   me.correct = "races")
                    ) %>%
    as_tibble() %>%
    arrange(id)


extr_probs = function(d) {
    arrange(d, id) %>%
        select(starts_with("pred_")) %>%
        `colnames<-`(c("pr_white", "pr_black", "pr_hisp", "pr_asian", "pr_other"))
}

r_probs = predict_race_sgz(last_name, zip, data=d_fit,
                           p_r=p_r, iterate=0, regularize=T)

fit0 = extr_probs(r_probs_naive) %>%
    model_race(party, zip, data=d_fit,
               reload_py=TRUE, config=list(n_mi=0, lr=0.25))
fit1 = extr_probs(r_probs_me) %>%
    model_race(party, zip, data=d_fit,
               reload_py=TRUE, config=list(n_mi=0, lr=0.25))
fit2 = model_race(r_probs$pr, party, zip, data=d_fit,
               reload_py=TRUE, config=list(n_mi=0, lr=0.25))

xr = list(
    true = p_xr,
    bisg = calc_joint_bisgz(extr_probs(r_probs_naive), d_fit$party),
    bisg_me = calc_joint_bisgz(extr_probs(r_probs_me), d_fit$party),
    bisgz = calc_joint_bisgz(r_probs$pr, d_fit$party),
    fit0 = calc_joint_model(fit0$p_xr, 0.5, p_r),
    fit0_low = calc_joint_model(fit0$p_xr, 0.05, p_r),
    fit0_high = calc_joint_model(fit0$p_xr, 0.95, p_r),
    fit1 = calc_joint_model(fit1$p_xr, 0.5, p_r),
    fit1_low = calc_joint_model(fit1$p_xr, 0.05, p_r),
    fit1_high = calc_joint_model(fit1$p_xr, 0.95, p_r),
    fit2 = calc_joint_model(fit2$p_xr, 0.5, p_r),
    fit2_low = calc_joint_model(fit2$p_xr, 0.05, p_r),
    fit2_high = calc_joint_model(fit2$p_xr, 0.95, p_r)
)

eval_joints(xr$true, "tv",
            bisg_cty_naive = xr$bisg,
            bisg_cty_me = xr$bisg_me,
            bisg_zip_naive = xr$bisgz,
            fit_cty_naive = xr$fit0,
            fit_cty_me = xr$fit1,
            fit_zip_naive = xr$fit2)

with(xr, mean(fit0_low <= true & true <= fit0_high))
with(xr, mean(fit1_low <= true & true <= fit1_high))
with(xr, mean(fit2_low <= true & true <= fit2_high))

print_cond(xr$true)
print_cond(xr$fit0 - xr$true)
print_cond(xr$fit1 - xr$true)
print_cond(xr$fit2 - xr$true)
print_cond(abs(xr$fit0 - xr$true) - abs(xr$fit1 - xr$true))

