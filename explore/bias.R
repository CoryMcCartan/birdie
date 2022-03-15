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
