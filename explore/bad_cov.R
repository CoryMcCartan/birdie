# library(raceproxy)
library(dplyr)
library(here)
library(ggplot2)
library(wacolors)

if (file.exists(voterfile <- here("data-raw/nc_voters.rds"))) {
    voters = readRDS(voterfile)
} else {
    voters = make_nc_statewide(voterfile)
}

d_fit = slice_sample(voters, n=20e3) %>%
    mutate(gender = as.factor(coalesce(gender, "U"))) %>%
    filter(!is.na(age))

p_xr = prop.table(table(d_fit$party, d_fit$race))
p_x = rowSums(p_xr)
p_r = colSums(p_xr)

bisg = list(
    no = predict_race_sgz(last_name, zip, data=d_fit, p_r=p_r, iterate=0),
    me = predict_race_sgz_me(last_name, zip, data=d_fit, p_r=p_r)
)

config=list(it_avg=400, tol_rhat=1.1)
fits = list(
    no_0  = model_race(bisg$no, party, zip, data=d_fit, config=config),
    no_g  = model_race(bisg$no, party, zip, c(gender), data=d_fit, config=config),
    no_a  = model_race(bisg$no, party, zip, c(age), data=d_fit, config=config),
    no_ga = model_race(bisg$no, party, zip, c(gender, age), data=d_fit, config=config),
    me_0  = model_race(bisg$me, party, zip, data=d_fit, config=config),
    me_g  = model_race(bisg$me, party, zip, c(gender), data=d_fit, config=config),
    me_a  = model_race(bisg$me, party, zip, c(age), data=d_fit, config=config),
    me_ga = model_race(bisg$me, party, zip, c(gender, age), data=d_fit, config=config)
)

xr = list(
    true = p_xr,
    bisg_no = calc_joint_bisgz(bisg$no, d_fit$party),
    bisg_me = calc_joint_bisgz(bisg$me, d_fit$party)
)
xr = c(xr, lapply(fits, \(f) calc_joint_model(f$p_xr, 0.5, p_r)))

do.call(eval_joints, c(list(xr$true, "tv"), xr))

do.call(eval_joints, c(list(xr$true, "tv_col"), xr)) %>%
    tail(-3) %>%
    rename(TV=TV_COL) %>%
    tidyr::unnest_longer(TV) %>%
    tidyr::separate(method, into=c("correction", "covs"), sep="_") %>%
ggplot(aes(TV_id, TV, color=covs, shape=correction)) +
    geom_jitter(size=4) +
    scale_color_wa_d()

do.call(eval_joints, c(list(xr$true, "tv_row"), xr)) %>%
    tail(-3) %>%
    rename(TV=TV_ROW) %>%
    tidyr::unnest_longer(TV) %>%
    tidyr::separate(method, into=c("correction", "covs"), sep="_") %>%
ggplot(aes(TV_id, TV, color=covs, shape=correction)) +
    geom_jitter(size=4) +
    scale_color_wa_d()
