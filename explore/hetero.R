library(tidyverse)
library(collapse)

run = function(N = 100e3, b1 = 0.5, b0 = -0.3, xz0=0.5, xz1=0.5, p = 0.75) {
    z = rbinom(N, 1, p=p)
    x = if_else(z == 0, rbinom(N, 1, xz0), rbinom(N, 1, xz1))
    y = rbinom(N, 1, if_else(z == 0, 0.7 + b0*x, 0.3 + b1*x))

    m = summary(lm(y ~ x))

    joint = matrix(c(1-xz0, xz0, 1-xz1, xz1) * c(1-p, 1-p, p, p), nrow=2)
    cond = prop.table(joint, margin=1)
    x0_avg = cond[1, 2]*(0.3) + cond[1, 1]*(0.7)
    x1_avg = cond[2, 2]*(0.3 + b1) + cond[2, 1]*(0.7 + b0)

    qDF(list(b1=b1, b0=b0, p=p,
             x0_avg=x0_avg, x1_avg=x1_avg,
             est0 = m$coefficients[1, 1],
             est1 = m$coefficients[1, 1] + m$coefficients["x", 1]))
}


res = map_dfr(1:4000, function(i) {
    # run(1e3, runif(1, -0.3, 0.7), runif(1, -0.7, -0.3), runif(1))
    run(1e3, runif(1, -0.3, 0.7), runif(1, -0.7, -0.3),
        xz0=runif(1, 0.1, 0.9), xz1=runif(1, 0.1, 0.9), p=runif(1))
}) |> as_tibble()

ggplot(res, aes(x1_avg, est1-x1_avg, color=b0)) +
    geom_point() +
    geom_smooth(color="black")

ggplot(NULL, aes(x, y, color=as.factor(z), group=z)) +
    stat_summary(geom="point", fun.data=function(x) tibble(y=mean(x), size=length(x)/2000)) +
    geom_smooth(method=lm) +
    geom_smooth(aes(group=1, color=1), method=lm, color="black")

p = pt(res$t, df=10e3-2)
hist(p)
ks.test(p, "punif")
