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

d_fit = slice(voters, 1:10e3)

# tables
p_xr = prop.table(table(d_fit$party, d_fit$race))
p_x = rowSums(p_xr)
p_r = colSums(p_xr)

fit = model_race(party, last_name, zip, Z=c(age, gender), data=d_fit, p_r=p_r,
                 regularize=T, alpha=3,
                 #methods=c("bis", "bisg", "nonparam", "additive"),
                 methods=c("bis", "bisg"),
                 stan_method="vb",
                 iter=400, verbose=TRUE)

fit$pyro = est_additive_pyro(d$X, d$GZ_mat, d$GZ_var, d$pr_base,
                             iter=20e3, lr=0.25, subsamp=2048, tol=2.50,
                             reload=TRUE)

plot(fit$pyro$loss, type='l', ylim=c(min(c(1, fit$pyro$loss)), max(c(2, fit$pyro$loss))))
lines(seq_along(fit$pyro$r_hat)*100, fit$pyro$r_hat, col='red')

xr = calc_joints(p_xr, voters, fit)
eval_joints(xr$true, "tv",
            bis = xr$bis,
            bisg = xr$bisg,
            nonparam = xr$nonparam,
            additive = xr$additive,
            pyro = xr$pyro)
