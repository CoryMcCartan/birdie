suppressMessages({
    library(dplyr)
    devtools::load_all(".")
    library(here)
})

d = readRDS(here("data-raw/nc_voters_small.rds")) |>
    filter(!is.na(party))
p_r = with(d, prop.table(table(race)))

r_probs_vn = bisg(~ nm(last_name) + zip(zip), d, p_r=p_r)
r_probs_me = bisg_me(~ nm(last_name) + zip(zip), d, p_r=p_r, cores=4)
r_probs = r_probs_me

data = d |>
    mutate(zip = proc_zip(zip),
           lic = c("no license", "license")[lic + 1],
           n_voted=as.factor(n_voted)) |>
    left_join(census_race_geo_table("zcta", counts=FALSE), by=c("zip"="GEOID")) |>
    mutate(across(white:other, ~ coalesce(., p_r[cur_column()]))) |>
    select(party, zip, county, race, gender, age, n_voted, lic, white:other)

Y =  data$party

x0 = birdie(r_probs, Y ~ 1, data)
x1 = birdie(r_probs, Y ~ zip, data)
x2 = birdie(r_probs, Y ~ white + black + hisp + (1|zip), data, ctrl=birdie.ctrl(abstol=1e-5))

{
    par(mfrow=c(3,2), mar=c(2, 2, 0.2, 0.2))
    triangle = cbind(c(0, 0.5, 1), c(0, 1, 0))
    for (r in 1:6) {
        theta = t(coef(x2, subgroup=TRUE)[1:3, r, ])
        plot(rbind(triangle, triangle[1, ]), type="l", col="#777777")
        points(theta %*% triangle, cex=0.7, pch=16, col="#00000022")
        points(coef(x2)[1:3, r] %*% triangle, cex=2, col="red", pch=8)
        text(triangle, labels=toupper(rownames(coef(x1))[1:3]), cex=1.0)
        text(0.1, 0.9, labels=names(p_r[r]), cex=2)
    }
    par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
}



xr = list(
    true = with(d, prop.table(table(Y, race))),
    pool = coef(x0) %*% diag(colMeans(fitted(x0))),
    sat = coef(x1) %*% diag(colMeans(fitted(x1))),
    mmm = coef(x2) %*% diag(colMeans(fitted(x2))),
    # sat2 = coef(x1b) %*% diag(colMeans(fitted(x1b))),
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
colSums(abs(to_cond(xr$true) - to_cond(xr$mmm)))/2
colSums(abs(to_cond(xr$true) - to_cond(xr$ols)))/2
colSums(abs(to_cond(xr$true) - to_cond(xr$weight)))/2
colSums(abs(to_cond(xr$true) - to_cond(xr$thresh)))/2


# try out multistage regression
ctrl = birdie.ctrl(abstol=1e-4)
x3 = birdie(r_probs, lic ~ zip, data, ctrl=ctrl) |>
    fitted() |>
    birdie(party ~ lic + (1|zip), data, ctrl=ctrl)
