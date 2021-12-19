suppressMessages(library(here))
devtools::load_all(here("."))

if (file.exists(voterfile <- here("data/nc_voters.rds"))) {
    voters = readRDS(voterfile)
} else {
    voters = make_nc_statewide(voterfile)
}

d_fit = slice_sample(voters, n=50e3) %>%
    mutate(gender = as.factor(coalesce(gender, "U"))) %>%
    filter(!is.na(age))

# tables
p_xr = prop.table(table(d_fit$party, d_fit$race))
p_x = rowSums(p_xr)
p_r = colSums(p_xr)

r_probs = predict_race_sgz(last_name, zip, c(gender, age), data=d_fit,
                           p_r=p_r, iterate=0)

#fit0 = model_race(r_probs, party, zip, c(age, gender), data=d_fit,
#                  reload_py=TRUE, config=list(n_mi=0, lr=0.3))

fit0 = model_race(r_probs$pr, party, zip, c(gender, age), data=d_fit, sgz=r_probs,
                  reload_py=TRUE, config=list(n_mi=0, lr=0.25))
fit1 = model_race(r_probs$pr, party, zip, c(gender, age), data=d_fit, sgz=r_probs,
                  reload_py=TRUE, config=list(n_mi=7, lr=0.25))

plot(c(tail(fit1$loss, -100), fit1$loss2), type='l')

fit = fit1
xr = list(
    true = p_xr,
    bisgz = calc_joint_bisgz(r_probs$pr, d_fit$party),
    model = calc_joint_model(fit$p_xr, 0.5, p_r),
    model_low = calc_joint_model(fit$p_xr, 0.05, p_r),
    model_high = calc_joint_model(fit$p_xr, 0.95, p_r)
)
cat("Coverage =", mean((xr$model_low[,1:2] < p_xr[,1:2]) & (p_xr[,1:2] < xr$model_high[,1:2])), "\n")

eval_joints(xr$true, "tv",
            bisgz = xr$bisgz,
            model = xr$model)

print_cond(xr$true, "true")
print_cond(xr$bisgz, "BISGZ")
print_cond(xr$model, "median")
print_cond(xr$model_low, "low")
print_cond(xr$model_high, "high")
print_cond(xr$model_high - xr$model_low, "difference")
xr$model_low < xr$true & xr$true < xr$model_high
