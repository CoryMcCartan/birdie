suppressMessages({
    library(dplyr)
    devtools::load_all(".")
    library(here)
})

d = readRDS(here("data-raw/nc_voters_small.rds")) |>
    filter(!is.na(party))
p_r = with(d, prop.table(table(race)))

r_probs = predict_race_sgz(last_name, zip, data=d, p_r=p_r,
                           est_r_gz=TRUE, iterate=0)
p_r_est = colMeans(r_probs)
alpha = c(10, 10, 10, 1)

# x = birdie(r_probs, party ~ 1, d, alpha=alpha)
data = mutate(d, zip = coalesce(zip, "<none>")) |>
    select(party, zip, county, race)
# x = birdie(r_probs, party ~ (1 | zip), data, alpha=alpha, iter=200)
x = birdie(r_probs, party ~ (1 | zip), data, alpha=alpha, iter=20)

xr = list(
    true = with(d, prop.table(table(party, race))),
    weight = calc_joint_bisgz(r_probs, d$party, "weight"),
    ols = calc_joint_bisgz(r_probs, d$party, "ols"),
    post_weight = calc_joint_bisgz(x$p_ryxs, d$party, "weight"),
    # map0 = x$map0 %*% diag(p_r_est),
    EM = x$map %*% diag(colMeans(x$p_ryxs)),
    EM2 = x$map %*% diag(p_r)
)
# xr = c(xr[1], lapply(xr[-1], \(tbl) rake(tbl, rowSums(xr$true), colSums(xr$true))))
xr = c(xr[1], lapply(xr[-1], \(tbl) tbl %*% diag(p_r / colSums(tbl))))

do.call(eval_joints, c(list(xr$true, "tv"), xr))

print_cond(xr$true)
print_cond(xr$EM)
print_cond(xr$ols)

colSums(abs(to_cond(xr$true) - to_cond(xr$EM)))/2
colSums(abs(to_cond(xr$true) - to_cond(xr$ols)))/2
colSums(abs(to_cond(xr$true) - to_cond(xr$weight)))/2
