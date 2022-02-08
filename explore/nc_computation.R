suppressMessages(library(here))
devtools::load_all(here("."))

if (file.exists(voterfile <- here("data/nc_voters.rds"))) {
    voters = readRDS(voterfile)
} else {
    counties = tigris::fips_codes$county[tigris::fips_codes$state == "NC"] %>%
        stringr::str_sub(end=-8)
    voters = do.call(rbind, lapply(counties, make_nc_df))
    voters = voters %>%
        select(last_name, party, race, zip, gender, age, birth_state, lic) %>%
        mutate(across(where(is.character), factor))
    saveRDS(voters, voterfile, compress="xz")
}

d_fit = slice_sample(voters, n=500e3) %>%
    mutate(gender = as.factor(coalesce(gender, "U"))) %>%
    filter(!is.na(age))

# tables
p_xr = prop.table(table(d_fit$party, d_fit$race))
p_x = rowSums(p_xr)
p_r = colSums(p_xr)

fit = model_race(party, last_name, zip, Z=c(age, gender), data=d_fit, p_r=p_r,
                 regularize=T, alpha=3,
                 #methods=c("bis", "bisg", "nonparam", "additive"),
                 methods=c("bis", "bisg", "pyro"),
                 stan_method="hmc",
                 iter=100, verbose=TRUE)


noise_pr = function(x, amt=0.2) {
    x = exp(log(x) + rnorm(length(x), 0, amt))
    x / rowSums(x)
}

for (i in 1:8) {
    pr_base = noise_pr(d$pr_base)
    fit_pyro = est_additive_pyro(d$X, d$GZ_mat, d$GZ_var, pr_base,
                                 iter=20e3, lr=0.25, draws=200, subsamp=2048,  tol=1.35,
                                 reload=TRUE)
    if (i == 1) {
        fit$pyro = fit_pyro
        fit$pyro_pooled = fit_pyro
    } else {
        fit$pyro_pooled$p_xr = abind::abind(fit$pyro_pooled$p_xr, fit_pyro$p_xr, along=1)
    }
}

plot(fit$pyro$loss, type='l', ylim=c(min(c(1, fit$pyro$loss)), max(c(2, fit$pyro$loss))))
lines(seq_along(fit$pyro$r_hat)*100, fit$pyro$r_hat, col='red')

xr = calc_joints(p_xr, d_fit, fit)
eval_joints(xr$true, "tv",
            bis = xr$bis,
            bisg = xr$bisg,
            #nonparam = xr$nonparam,
            additive = xr$additive,
            pyro = xr$pyro)#,
            #pyro_pooled = xr$pyro_pooled)

print_cond(xr$true, "true")
print_cond(xr$pyro, "pyro")
print_cond(xr$pyro_low, "low")
print_cond(xr$pyro_high, "high")
