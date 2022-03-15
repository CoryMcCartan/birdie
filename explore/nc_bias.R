suppressMessages(library(here))
devtools::load_all(here("."))

if (file.exists(voterfile <- here("data/nc_voters.rds"))) {
    voters = readRDS(voterfile)
} else {
    voters = make_nc_statewide(voterfile)
}

d_fit = slice_sample(voters, n=10e3) %>%
    mutate(gender = as.factor(coalesce(gender, "U"))) %>%
    filter(!is.na(age))

# tables
p_xr = prop.table(table(d_fit$party, d_fit$race))
p_x = rowSums(p_xr)
p_r = colSums(p_xr)

r_probs = predict_race_sgz(last_name, zip, c(gender, age), data=d_fit,
                           p_r=p_r, iterate=0)

fit = model_race(r_probs$pr, party, zip, c(gender, age), data=d_fit, sgz=r_probs,
                  reload_py=TRUE, config=list(n_mi=0, lr=0.25))

xr = list(
    true = p_xr,
    bisgz = calc_joint_bisgz(r_probs$pr, d_fit$party),
    model = calc_joint_model(fit$p_xr, 0.5, p_r),
    model_low = calc_joint_model(fit$p_xr, 0.05, p_r),
    model_high = calc_joint_model(fit$p_xr, 0.95, p_r)
)
cat("Coverage =", mean((xr$model_low < p_xr) & (p_xr < xr$model_high)), "\n")

eval_joints(xr$true, "tv",
            bisgz = xr$bisgz,
            model = xr$model)
