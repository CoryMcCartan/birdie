if (file.exists(voterfile <- here("data-raw/nc_voters.rds"))) {
    voters = readRDS(voterfile)
} else {
    voters = make_nc_statewide(voterfile)
}

d = slice_sample(voters, n=200e3) %>%
    mutate(gender = as.factor(coalesce(gender, "U"))) %>%
    filter(!is.na(age))

p_xr = prop.table(table(d$party, d$race))
p_x = rowSums(p_xr)
p_r = colSums(p_xr)

# read in census
d_cens = readRDS(system.file("extdata", "zip_race_2020_nc.rds",
                            package="raceproxy", mustWork=TRUE)) %>%
    select(-starts_with("pop"))
for (i in seq_along(p_r)) {
    d_cens[, 2+i] = d_cens[, 2+i] / d_cens[, 2]
}
#d_cens = d_cens[, -2:-9]
d_cens = d_cens[, -2]
colnames(d_cens) = c("zip", "white", "black", "hisp", "asian", "aian", "other")

bisg_10 = predict_race_sgz(last_name, zip, data=d, p_r=p_r)
bisg_20 = predict_race_sgz(last_name, zip, data=d, p_r=p_r, p_rgz=d_cens)


fit_10 = model_race(bisg_10, party, zip, gender, data=d, config=list(tol_rhat=1.15))
fit_20 = model_race(bisg_20, party, zip, gender, data=d, config=list(tol_rhat=1.15))

xr = list(
    true = p_xr,
    bisg_10 = calc_joint_bisgz(bisg_10, d$party),
    bisg_20 = calc_joint_bisgz(bisg_20, d$party),
    model_10 = calc_joint_model(fit_10$p_xr, 0.5, p_r),
    model_20 = calc_joint_model(fit_20$p_xr, 0.5, p_r)
)

do.call(eval_joints, c(list(xr$true, "tv"), xr))
do.call(eval_joints, c(list(xr$true, "tv_col"), xr)) %>%
    tidyr::unnest_longer(TV_COL) %>%
    pivot_wider(names_from=TV_COL_id, values_from=TV_COL)
do.call(eval_joints, c(list(xr$true, "tv_row"), xr)) %>%
    tidyr::unnest_longer(TV_ROW) %>%
    pivot_wider(names_from=TV_ROW_id, values_from=TV_ROW)
