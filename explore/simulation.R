library(tidyverse)
library(patchwork)
library(scales)
library(cli)

#' Generate data
#'
#' @param N observations
#' @param n_zx 2 integers: the number of Z (joint with race) and X (not joint) covariates
#' @param x_conf the amount of confounding: how well X predicts race
#' @param x_conf the amount of confounding: how well X predicts race
#' @param r_bias the amount of nonlinear distortion to apply to the race probabilities. Positive skews them upward.
#' @param y_rate the target marginal mean of the outcome Y
#'
#' @returns a tibble
gen_data = function(N=100, n_zx=c(3, 3), x_conf=0.0, z_pred=0.8, r_bias=0.0, y_rate=0.1) {
    K = sum(n_zx)
    ZX = matrix(as.integer(rbernoulli(N * K, 0.5)), nrow=N, ncol=K)
    colnames(ZX) = c(paste0("z_", 1:n_zx[1]), paste0("x_", 1:n_zx[2]))

    y_coef = c(seq(0.5, -0.5, length.out=n_zx[1]),
               seq(1, -1, length.out=n_zx[2]))
    r_coef = c(seq(1, -1, length.out=n_zx[1]),
               seq(1, -1, length.out=n_zx[2]) * x_conf)

    R_p = pnorm(z_pred * (ZX %*% r_coef) - 0.5)
    R = as.integer(rbernoulli(N, R_p))
    R_est = fitted(glm.fit(ZX[, 1:n_zx[1]], R, family=binomial("probit"))) ^ exp(-r_bias)

    Y_linpred = ZX %*% y_coef
    Y_p = pnorm(Y_linpred - mean(Y_linpred) + qnorm(y_rate))
    Y = as.integer(rbernoulli(N, Y_p))
    Y_est = fitted(glm.fit(cbind(1, ZX), Y, family=binomial("probit")))

    cbind(Y=Y, R=R, ZX, Y_p=Y_p[, 1], Y_est=Y_est, R_p=R_p[, 1], R_est=R_est) |>
        as_tibble()
}

#' Plots estimated vs true probabilities for Y and R
#'
#' @inheritParams gen_data
ex_plot = function(N=1000, n_zx=c(5, 10), ...) {
    d = gen_data(N, n_zx, ...)

    p1 = ggplot(d, aes(Y_p, Y_est)) +
        geom_count(alpha=0.2) +
        geom_abline(slope=1, color="red") +
        coord_equal(xlim=0:1, ylim=0:1) +
        scale_x_continuous("Pr(Y)", labels=percent) +
        scale_y_continuous("Est. Y", labels=percent) +
        labs(title="Outcome (Y)") +
        guides(size="none")
    p2 = ggplot(d, aes(R_p, R_est)) +
        geom_count(alpha=0.2) +
        geom_abline(slope=1, color="red") +
        coord_equal(xlim=0:1, ylim=0:1) +
        scale_x_continuous("Pr(R)", labels=percent) +
        scale_y_continuous("Est. R", labels=percent) +
        labs(title="Race (R)") +
        guides(size="none")
    p1 + p2
}

#' Calculate the imbalance between Y and R_est, given Y_est
#'
#' @param d from `gen_data()`
#' @param n_bin the number of bins for Y_est
#'
#' @return a numeric
cond_dev = function(d, n_bin=5) {
    bins = cut(d$Y_est, n_bin)
    mean(abs(tapply(d$R_est[d$Y == 1], bins[d$Y == 1], mean) -
                 tapply(d$R_est[d$Y == 0], bins[d$Y == 0], mean)))
}

#' Run a simulation
#'
#' @param N_sim the number of simulated datasets
#' @inheritParams gen_data
#'
#' @returns a tibble with the following entries for each simulated dataset:
#' * `diff_act` the true racial disparity: E[Y|R=1] - E[Y|R=0]
#' * `diff_p` the racial disparity calculated using the true R probabilities
#' * `diff_est` the racial disparity calculated using the estimated R probabilities
#' * `indep_p` the p value for R_est in a GLM of Y on R_est and Y_est
#' * `dev` the imbalance metric calculatec from `cond_dev()`
run_sim = function(N_sim=500, N=1000, n_zx=c(4, 8), ...) {
    map_dfr(cli_progress_along(1:N_sim), function(i) {
        d = gen_data(N, n_zx, ...)
        m = glm(Y ~ R_est + Y_est, data=d, family=binomial("probit"))

        fn = function(b) {
            m = glm(Y ~ I(R_est^exp(b)) + Y_est, data=d,
                    family=binomial("probit"))
            abs(summary(m)$coefficients[2, 1])
            #abs(summary(m)$coefficients[2, 4] - 0.5)
        }
        opt_res = optimize(fn, c(-2, 2))
        b = exp(opt_res$minimum)
        m_opt = glm(Y ~ I(R_est^b) + Y_est, data=d,
                    family=binomial("probit"))

        tibble(diff_act = diff(tapply(d$Y, d$R, mean)),
               diff_p = weighted.mean(d$Y, d$R_p) - weighted.mean(d$Y, 1 - d$R_p),
               diff_est = weighted.mean(d$Y, d$R_est) - weighted.mean(d$Y, 1 - d$R_est),
               diff_corr = weighted.mean(d$Y, d$R_est^b) - weighted.mean(d$Y, 1 - d$R_est^b),
               indep_p = summary(m)$coefficients["R_est", 4],
               corr_p = summary(m_opt)$coefficients[2, 4],
               corr_b = b,
               dev = cond_dev(d))
    })
}

#' RMSE
rmse = function(x, y, data=NULL) {
    x = rlang::eval_tidy(rlang::enquo(x), data)
    y = rlang::eval_tidy(rlang::enquo(y), data)
    sqrt(mean((x - y)^2))
}
