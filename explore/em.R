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
r_probs = r_probs_me

p_r_est = colMeans(r_probs)

data = mutate(d, zip = proc_zip(zip), n_voted=as.factor(n_voted)) |>
    select(party, zip, county, race, gender, age, n_voted)

if (FALSE) {
    formula = party ~ zip
    ctrl = birdie.ctrl()
    prior = rep(1.0001, 4)
    p_rxs = as.matrix(r_probs)
    # later
    Y = Y_vec

    ests = dirichlet_map(as.integer(d$party), rep_len(1, nrow(d)), as.matrix(r_probs), rep(1.001, 4), 1)
    em_step = function(curr) em_dirichlet(curr, as.integer(d$party), rep_len(1, nrow(d)), as.matrix(r_probs), rep(1.001, 4), 1)
}

Y =  data$n_voted

x0 = birdie(r_probs, Y ~ 1, data)
x1 = birdie(r_probs, Y ~ zip, data)
# x = birdie(r_probs, party ~ (1 | zip), data, ctrl=birdie.ctrl(max_iter=20))

xr = list(
    true = with(d, prop.table(table(Y, race))),
    pool = x0$map %*% diag(colMeans(x0$p_ryxs)),
    sat = x1$map %*% diag(colMeans(x1$p_ryxs)),
    # glmm = x$map %*% diag(colMeans(x$p_ryxs)),
    # pols = calc_joint_bisgz_ols(r_probs, d$party, d$zip, with(d, prop.table(table(zip, race), 2))),
    ols = calc_joint_bisgz(r_probs, Y, "ols"),
    weight = calc_joint_bisgz(r_probs, Y, "weight"),
    thresh = calc_joint_bisgz(r_probs, Y, "thresh")
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


# try out multistage regression
x0 = birdie(r_probs, age ~ zip, data)
x1 = birdie(x0$p_ryxs, party ~ zip * age, data)
x1 = birdie(x0$p_ryxs, party ~ age, data)
