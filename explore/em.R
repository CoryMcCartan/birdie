suppressMessages({
    library(dplyr)
    devtools::load_all(".")
    library(here)
})

d = readRDS(here("data-raw/nc_voters_small.rds")) |>
    filter(!is.na(party))
p_r = with(d, prop.table(table(race)))

r_probs = bisg(~ nm(last_name) + zip(zip), d, p_r=p_r)
r_probs_me = bisg_me(~ nm(last_name) + zip(zip), d, p_r=p_r, cores=4)

p_r_est = colMeans(r_probs)

data = mutate(d, zip = proc_zip(zip)) |>
    select(party, zip, county, race, n_voted)

if (FALSE) {
    formula = party ~ zip
    ctrl = birdie.ctrl()
    prior = rep(1.0001, 4)
    p_rxs = as.matrix(r_probs)
    # later
    Y = Y_vec
}

x0 = birdie(r_probs, party ~ 1, data, prior=rep(1.01, 4))
x1 = birdie(r_probs, party ~ zip, data, prior=rep(1.01, 4))
# x = birdie(r_probs, party ~ (1 | zip), data, ctrl=birdie.ctrl(max_iter=20))

xr = list(
    true = with(d, prop.table(table(party, race))),
    pool = x0$map %*% diag(colMeans(x0$p_ryxs)),
    sat = x1$map %*% diag(colMeans(x1$p_ryxs)),
    # glmm = x$map %*% diag(colMeans(x$p_ryxs)),
    # pols = calc_joint_bisgz_ols(r_probs, d$party, d$zip, with(d, prop.table(table(zip, race), 2))),
    ols = calc_joint_bisgz(r_probs, d$party, "ols"),
    weight = calc_joint_bisgz(r_probs, d$party, "weight"),
    thresh = calc_joint_bisgz(r_probs, d$party, "thresh")
)

# xr = c(xr[1], lapply(xr[-1], \(tbl) rake(tbl, rowSums(xr$true), colSums(xr$true))))
xr = c(xr[1], lapply(xr[-1], \(tbl) tbl %*% diag(p_r / colSums(tbl))))

do.call(eval_joints, c(list(xr$true, "tv"), xr))

print_cond(xr$true)
print_cond(xr$pool)
print_cond(xr$sat)
print_cond(xr$ols)

colSums(abs(to_cond(xr$true) - to_cond(xr$pool)))/2
colSums(abs(to_cond(xr$true) - to_cond(xr$sat)))/2
# colSums(abs(to_cond(xr$true) - to_cond(xr$glmm)))/2
colSums(abs(to_cond(xr$true) - to_cond(xr$ols)))/2
colSums(abs(to_cond(xr$true) - to_cond(xr$weight)))/2
