suppressMessages({
    library(dplyr)
    devtools::load_all(".")
    library(here)
})

d = readRDS(here("data-raw/nc_voters_small.rds")) |>
    filter(!is.na(party)) |>
    head(1e4)
p_r = with(d, prop.table(table(race)))
p_xr = with(d, prop.table(table(party, race), 2))

r_probs_vn = bisg(~ nm(last_name) + zip(zip), d, p_r=p_r)
# r_probs_me = bisg_me(~ nm(last_name) + zip(zip), d, p_r=p_r, cores=4)
r_probs = r_probs_vn

# do birdie -------

data = d |>
    mutate(zip = proc_zip(zip),
           lic = c("no license", "license")[lic + 1],
           n_voted=as.factor(n_voted)) |>
    left_join(census_race_geo_table("zcta", counts=FALSE), by=c("zip"="GEOID")) |>
    mutate(across(white:other, ~ coalesce(., p_r[cur_column()]))) |>
    select(party, zip, county, race, gender, age, n_voted, lic, white:other)
data$Yc = (data$party == "dem") + rnorm(nrow(data))

Y =  data$party

xw = est_weighted(r_probs, Y ~ 1, data)
x0 = birdie(r_probs, Y ~ 1, data)
x1 = birdie(r_probs, Y ~ zip, data, algorithm="em_boot")
x1b = birdie(r_probs, Y ~ zip, data, algorithm="gibbs", iter=1e3)
x2 = birdie(r_probs, Y ~ white + black + hisp + (1|zip), data,
            prior = list(scale_int=2, scale_sigma=0.05, scale_beta=0.2),
            family = cat_mixed(), ctrl=birdie.ctrl(abstol=5e-6))
x3 = birdie(r_probs, Yc ~ black + white, data, family = gaussian(), ctrl=birdie.ctrl(abstol=5e-6))
x3b = birdie(r_probs, Yc ~ black + white, data, family = gaussian(), ctrl=birdie.ctrl(abstol=5e-6))
x2 = birdie(r_probs, party=='dem' ~ black + white, data, family = gaussian(),
            ctrl=birdie.ctrl(abstol=5e-6),
            prior=list(scale_beta=2.5, scale_int=5, n_sigma=20000, loc_sigma=0.5))
x3 = birdie(r_probs, party == "dem" ~ 1, data, family=gaussian(),
            algorithm="gibbs", ctrl=birdie.ctrl(abstol=1e-4, fix_sigma=TRUE))


{
    par(mfrow=c(3,2), mar=c(2, 2, 0.2, 0.2))
    triangle = cbind(c(0, 0.5, 1), c(0, 1, 0))
    for (r in 1:6) {
        theta = t(coef(x2, subgroup=TRUE)[1:3, r, ])
        plot(rbind(triangle, triangle[1, ]), type="l", col="#777777")
        points(theta %*% triangle, cex=0.7, pch=16, col="#00000022")
        points(coef(x2)[1:3, r] %*% triangle, cex=2, col="red", pch=8)
        text(triangle, labels=toupper(rownames(coef(x1))[1:3]), cex=1.0)
        text(0.1, 0.9, labels=toupper(names(p_r[r])), cex=1.5)
    }
    par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
}

sum(abs(coef(x1) - p_xr) %*% p_r) / 2
sum(abs(coef(x1b) - p_xr) %*% p_r) / 2

sum(colSums(abs(coef(x0) - p_xr)) * p_r) / 2
sum(colSums(abs(coef(x1) - p_xr)) * p_r) / 2
sum(colSums(abs(coef(x2) - p_xr)) * p_r) / 2


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

# try gibbs --------
Y = as.integer(data$party)
X = rep_along(Y, 1)
R0 = as.integer(simulate(r_probs))
prior = check_make_prior_cat_dir(NULL, data$party, p_rxs, levels(data$race))

N_sim = 500
ests = matrix(nrow=4*6, ncol=N_sim)

p_rxs = as.matrix(r_probs)
ests[, 1] = c(prop.table(1 + table(R0, Y), 1))
# ests = dirichlet_map(Y, X, p_rxs, prior$alpha, 1)

for (i in 2:N_sim) {
    p_ryxs = calc_bayes(Y, X, ests[, i-1], p_rxs, 1, 4)
    R = mat_rcatp(p_ryxs)
    alpha = 1 + table(R, Y)
    r_ests = rgamma(length(alpha), alpha, 1) |>
        matrix(nrow=6, ncol=4)
    ests[, i] = c(r_ests / c(r_ests %*% rep(1, 4)))
}

matplot(t(ests[6*0 + 1:6, ]), type='l', lty="solid",
        col=c("black", gray.colors(6)[1:5]), xlab="Iteration", ylab="Pr(party | white)")
matrixStats::rowSds(ests) |> matrix(nrow=4, ncol=6, byrow=TRUE)
ests_final = matrix(rowMeans(ests), nrow=4, ncol=6, byrow=TRUE)
ests_em = coef(birdie(r_probs, party ~ 1, data))
round(abs(ests_em - ests_final), 2)

sum(abs(ests_em - p_xr) %*% p_r) / 2
sum(abs(ests_final - p_xr) %*% p_r) / 2

vcov = chol2inv(qr.R(res$qr)) * sum(res$residuals^2) / (nrow(p_rxs) - ncol(p_rxs))
sigma_draw = sum(res$residuals^2) / rchisq(1, df=nrow(p_rxs) - ncol(p_rxs))
beta_draw = res$coefficients + backsolve(qr.R(res$qr), diag(6)) %*% rnorm(6, sd=sigma_draw)
bench::mark(
    rowSums(r_ests),
    as.vector(r_ests %*% rep(1, 4)),
    c(r_ests %*% rep(1, 4)),
)


library(dbarts)
Y = 1*(data$party == "dem")
mb = dbarts(Yc ~ R0 + gender + white + black, data=data)


# try out multistage regression -----
ctrl = birdie.ctrl(abstol=1e-4)
x3 = birdie(r_probs, lic ~ zip, data, ctrl=ctrl) |>
    fitted() |>
    birdie(party ~ lic + (1|zip), data, ctrl=ctrl)


# marginalizing --------

d_gr = as_tibble(attr(r_probs, "p_rgx")) |>
    mutate(zip = x2$tbl_gx$zip) |>
    tidyr::pivot_longer(-zip, names_to="race", values_to="wt") |>
    group_by(zip) |>
    summarize(wt = sum(wt))

v1 = tidy(x2, subgroup=TRUE) |>
    left_join(d_gr, by=c("zip")) |>
    group_by(Y, race) |>
    summarize(estimate = weighted.mean(estimate, wt)) |>
    arrange(Y, race)

v2 = tidy(x2) |>
    arrange(Y, race)

x2b = est_weighted(fitted(x2), Y ~ (zip == "<none>"), data)

v3 = tidy(x2b) |>
    arrange(Y, race)

all.equal(v2, v3)


# OLS - birdie connection ------
X = as.matrix(r_probs)
y = as.numeric(d$party == "dem")
U = solve(t(X) %*% X)
W = nrow(X) * X %*% U # implied ols weighting matrix
idx = sample(nrow(X), 1e3)

Xb = birdie(r_probs, as.factor(party == "dem") ~ 1, data=d) |>
    fitted() |>
    as.matrix()

cbind(X[idx, 1], Xb[idx, 1], W[idx, 1]) |>
    pairs(labels=c("BISG prob.", "BIRDiE prob.", "OLS weight"), main="White")
cbind(X[idx, 2], Xb[idx, 2], W[idx, 2]) |>
    pairs(labels=c("BISG prob.", "BIRDiE prob.", "OLS weight"), main="Black")
cbind(X[idx, 3], Xb[idx, 3], W[idx, 3]) |>
    pairs(labels=c("BISG prob.", "BIRDiE prob.", "OLS weight"), main="Hispanic")

cbind(X[idx, 5], Xb[idx, 5], W[idx, 5]) |>
    cor()

plot(X[idx, 1], W[idx, 1])
plot(X[idx, 1], Xb[idx, 1])
plot(X[idx, 2], W[idx, 2])
weighted.mean(y, W[, 1])
weighted.mean(y, X[, 1])
lm.fit(X, y)$coefficients

colSums(y %*% W / nrow(X))
colSums(y %*% Xb) / colSums(Xb)

pairs(X[idx, ], cex=0.1)
pairs(W[idx, ], cex=0.1)
pairs(t(beta1)[idx, ], cex=0.1)

est = lm.fit(X, y)$coefficients #* exp(rnorm(6)*0.5)
names(est) = levels(d$race)
# est = 0*est + 0.5
y2 = rnorm(nrow(X), est[d$race], 0.1)
# y2 = rnorm(nrow(X))
sigma = sd(y2 - X %*% est)
v_ls = -0.5 * ((y2 - X %*% est) / sigma)^2
v_br = log(rowSums(exp(-0.5 * (outer(y2, est, `-`) / sigma)^2) * X))
plot(v_br[idx], v_ls[idx], cex=0.3, xlab="log-lik BIRDiE", ylab="log-lik OLS")
cor(v_br, v_ls)
plot(y2[idx], v_ls[idx], cex=0.2, xlab="obs", ylab="log-lik")
points(y2[idx], v_br[idx], cex=0.1, col="red")

d |>
    mutate(.resid = residuals(lm(W[, 1] ~ X[, 1]))) |>
    arrange(desc(abs(.resid)))


theta = c(-.2, .2)
p = c(0.8, 0.2)
curve(-0.5*(x - sum(p * theta))^2, -5, 5)
curve(Vectorize(\(x) log(sum(p * exp(-0.5*(x - theta)^2))))(x), -5, 5, add=T, col='red')
