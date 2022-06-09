library(dplyr)
library(here)
devtools::load_all(here())
# library(ggplot2)
# library(wacolors)

# if (file.exists(voterfile <- here("data-raw/nc_voters.rds"))) {
#     voters = readRDS(voterfile)
# } else {
#     voters = make_nc_statewide(voterfile)
# }
#
# d_fit = slice_sample(voters, n=20e3) %>%
#     mutate(gender = as.factor(coalesce(gender, "U"))) %>%
#     filter(!is.na(age))
d_fit = readRDS("~/Desktop/voters.rds")

p_xr = prop.table(table(d_fit$party, d_fit$race))
p_x = rowSums(p_xr)
p_r = colSums(p_xr)

# bisg = predict_race_sgz_me(last_name, zip, data=d_fit, p_r=p_r)
bisg = predict_race_sgz(last_name, zip, data=d_fit, p_r=p_r)
# fit = model_race(bisg, party, zip, Z=c(gender, age), condition=gender,
#                  data=d_fit, config=list(it_avg=400, tol_rhat=1.1))
fit = model_race(bisg, party, zip, Z=c(gender, age), condition=gender,
                 data=d_fit, config=list(it_avg=100, tol_rhat=1.3))


xr = list(
    true = p_xr,
    weight = calc_joint_bisgz(bisg, d_fit$party, method="weight"),
    thresh = calc_joint_bisgz(bisg, d_fit$party, method="thresh"),
    mi = calc_joint_bisgz(bisg, d_fit$party, method="mi"),
    ols = calc_joint_bisgz(bisg, d_fit$party, method="ols"),
    model = calc_joint_model(fit$draws$global, 0.5, p_r)
)

xr_f = list(
    true = prop.table(table(d_fit$party[d_fit$gender=="F"], d_fit$race[d_fit$gender=="F"])),
    weight = calc_joint_bisgz(bisg[d_fit$gender == "F", ], d_fit$party[d_fit$gender == "F"], method="weight"),
    model_f = calc_joint_model(fit$draws$`gender: F`, 0.5, p_r)
)
xr_m = list(
    true = prop.table(table(d_fit$party[d_fit$gender=="M"], d_fit$race[d_fit$gender=="M"])),
    weight = calc_joint_bisgz(bisg[d_fit$gender == "M", ], d_fit$party[d_fit$gender == "M"], method="weight"),
    model_m = calc_joint_model(fit$draws$`gender: M`, 0.5, p_r)
)

do.call(eval_joints, c(list(xr$true, "tv"), xr))
do.call(eval_joints, c(list(xr_f$true, "tv"), xr_f))
do.call(eval_joints, c(list(xr_m$true, "tv"), xr_m))
