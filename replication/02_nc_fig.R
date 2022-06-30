fits_party <- read_rds(here("data-out/nc_fits_party.rds"))
fits_turnout <- read_rds(here("data-out/nc_fits_turnout.rds"))

xr = list(true = p_xr)
xr = c(xr, flatten(imap(fits, function(f, level) {
    out = list()
    out[[str_c("model_", level)]] = calc_joint_model(f, "global", 0.5, p_r)
    out
})))
xr = c(xr, flatten(imap(r_probs, function(d_pr, level) {
    out = list()
    out[[str_c("weight_", level)]] = calc_joint_bisgz(d_pr, d$party, method="weight")
    out[[str_c("thresh_", level)]] = calc_joint_bisgz(d_pr, d$party, method="thresh")
    out[[str_c("ols_", level)]] = calc_joint_bisgz(d_pr, d$party, method="ols")
    out
})))

do.call(eval_joints, c(list(xr$true, "tv"), xr))

do.call(eval_joints, c(list(xr$true, "tv_col"), xr)) %>%
    rename(TV=TV_COL) %>%
    unnest_longer(TV) %>%
    pivot_wider(names_from=TV_id, values_from=TV) %>%
    arrange(white)
