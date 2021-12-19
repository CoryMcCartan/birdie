library(dplyr)
library(here)
library(ggplot2)
library(wacolors)

if (file.exists(voterfile <- here("data-raw/nc_voters.rds"))) {
    voters = readRDS(voterfile)
} else {
    voters = make_nc_statewide(voterfile)
}

d_fit = slice_sample(voters, n=50e3) %>%
    mutate(gender = as.factor(coalesce(gender, "U"))) %>%
    filter(!is.na(age))

p_xr = prop.table(table(d_fit$party, d_fit$race))
p_x = rowSums(p_xr)
p_r = colSums(p_xr)

bisg = predict_race_sgz_me(last_name, zip, data=d_fit, p_r=p_r)
fit = model_race(bisg, party, zip, data=d_fit, config=list(it_avg=400, tol_rhat=1.1))


xr = list(
    true = p_xr,
    weight = calc_joint_bisgz(bisg, d_fit$party, method="weight"),
    thresh = calc_joint_bisgz(bisg, d_fit$party, method="thresh"),
    mi = calc_joint_bisgz(bisg, d_fit$party, method="mi"),
    ols = calc_joint_bisgz(bisg, d_fit$party, method="ols"),
    model = calc_joint_model(fit$p_xr, 0.5, p_r)
)

do.call(eval_joints, c(list(xr$true, "tv"), xr))
