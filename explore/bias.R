library(gtools)
library(mvtnorm)
library(tidyverse)

n_x = 4
n_r = 5
p_xr = matrix(rdirichlet(1, rep(2, n_x*n_r)), nrow=n_x, ncol=n_r)
p_xr = 0.05 + 0.9 * p_xr %*% diag(1 / colSums(p_xr))
p_r = rdirichlet(1, rep(1, n_r))[1, ]

theta = rmvnorm(200, as.numeric(p_xr), diag(0.05^2, length(p_xr)))
dim(theta) = c(200, n_x, n_r)

bound = function(N = 10, eps = 0.01) {
    x = sample(n_x, N, replace=T, prob=rowMeans(p_xr))
    r = sample(n_r, N, replace=T, prob=p_r)
    rhat = matrix(rdirichlet(N, 10*p_r), nrow=N, ncol=n_r)

    theta_i = theta[, x, ]
    d1 = sapply(1:N, \(i) theta_i[, i, ] %*% rhat[i, ])
    d2 = apply(theta_i, c(1, 2), \(x) sqrt(sum(x^2)))

    inner = sapply(1:N, \(i) {
        f2 = colMeans(theta_i[, i, ] / d2[, i])
        rowSums(sweep(theta_i[, i, ] / d1[, i], 2, f2)^2)
    })

    sqrt(mean(expm1(rowSums(inner)*eps^2))) / sqrt(N)
}

rho = 0.90
{
x = rmvnorm(1e6, c(0.5, 0.5), 0.05^2 * (rho + (1 - rho)*diag(2)))
#cor(x[,1], x[,2], method="pearson")
cor(log(x[,1]), x[,2], method="pearson")
cor(log(x[,1]), log(x[,2]), method="pearson")
cor(x[,1], x[,2], method="spearman")
}

library(mvtnorm)
library(gtools)
fn = function(...) {
    x = rmvnorm(1, sigma=diag(10))[1, ]
    x = (x - mean(x)) / sum(x^2)
    r = rdirichlet(1, rep(1, 10))
    x = x * r

    if (any(r - x < 0)) stop()

    out = r * sum(x^2) - x * sum(r * x)
    all(out * sign(sum(r*x)) >= 0)
}
res = map_dbl(1:5000, fn)
hist(res, breaks=100)
hist(sign(res), breaks=100)




# OLS and distribution of X ----

fn = function(loc=0.5, ...) {
    N = 100
    n_beta_obs = 50
    X = rbeta(N, loc*n_beta_obs, (1 - loc)*n_beta_obs)
    pr = 0.4*(1-X) + 0.8*X
    Y = rbinom(N, 1, pr)
    coef = lm.fit(cbind(1-X, X), Y)$coefficients
    tibble(loc = loc,
           beta1 = coef[1],
           beta2 = coef[2],
           mean_error = mean(Y - pr))
}

res = seq(0.01, 0.99, length.out=99) %>%
    sample(5000, replace=T) %>%
    map_dfr(fn)

ggplot(res, aes(loc, beta1)) +
    geom_point() +
    geom_smooth(formula=y~s(x), method="gam")

