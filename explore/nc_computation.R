suppressMessages(library(here))
devtools::load_all(here("."))

if (file.exists(voterfile <- here("data/nc_voters.rds"))) {
    voters = readRDS(voterfile)
} else {
    voters = make_nc_statewide(voterfile)
}

d_fit = slice_sample(voters, n=1e3) %>%
    mutate(gender = as.factor(coalesce(gender, "U"))) %>%
    filter(!is.na(age))

# tables
p_xr = prop.table(table(d_fit$party, d_fit$race))
p_x = rowSums(p_xr)
p_r = colSums(p_xr)

r_probs = predict_race_sgz(last_name, zip, Z=c(age, gender), data=d_fit,
                           p_r=p_r, iterate=0)

fit = model_race(r_probs, party, zip, c(age, gender), data=d_fit,
                 reload_py=TRUE)
plot(tail(fit$loss, -100), type='l')

xr = list(
    true = p_xr,
    bisgz = calc_joint_bisgz(r_probs, d_fit$party),
    model = calc_joint_model(fit$p_xr, 0.5, p_r),
    model_low = calc_joint_model(fit$p_xr, 0.05, p_r),
    model_high = calc_joint_model(fit$p_xr, 0.95, p_r)
)
cat("Coverage =", mean((p_xr > xr$model_low) & (p_xr < xr$model_high)), "\n")

eval_joints(xr$true, "tv",
            bisgz = xr$bisgz,
            model = xr$model)

print_cond(xr$true, "true")
print_cond(xr$bisgz, "BISGZ")
print_cond(xr$model, "median")
print_cond(xr$model_low, "low")
print_cond(xr$model_high, "high")
