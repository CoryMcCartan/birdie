library(tidyverse)
library(gtools)
library(cli)
library(scales)

gen_probs = function(n_r=2, n_x=5, n_s=3, s_strength) {
    P_r = as.numeric(rdirichlet(1, rep(2, n_r)))
    P_rs = rdirichlet(n_r, rep(1/s_strength, n_s))
    P_rs = diag(P_r) %*% P_rs
    P_xr = t(rdirichlet(n_r, rep(1, n_x)))
    list(r = P_r,
         rs = P_rs, # R, S joint
         xr = P_xr, # X | R
         xs = P_xr %*% P_rs)
}

samp_cond = function(N, x, cprob) {
    n_x = ncol(cprob)
    n_y = nrow(cprob)
    Y = integer(N)
    for (i in seq_len(n_x)) {
        idx = which(x == i)
        Y[idx] = sample.int(n_y, length(idx), replace=TRUE, prob=cprob[, i])
    }
    factor(Y, 1:n_y)
}


l_mad = function(x, y) mean(abs(x - y))
l_rmse = function(x, y) sqrt(mean((x - y)^2))

run_est = function(n_r=4L, n_x=100L, n_s=1e4L, N=1e3L, s_strength=1, loss=l_mad, impute=FALSE) {
    set.seed(5118)
    cli_h1("Confouding estimation with matrix inversion")
    cli_ul()
    cli_li("{comma(n_r)} race categories")
    cli_li("{comma(n_x)} X values")
    cli_li("{comma(n_s)} S values")
    cli_end()

    cli_h2("With true probabilities")
    probs = gen_probs(n_r, n_x, n_s, s_strength)
    P_r = probs$r
    P_rs = probs$rs
    P_xr = probs$xr
    P_xs = probs$xs

    P_xr_est = P_xs %*% MASS::ginv(P_rs)
    all.equal(P_xr_est, P_xr)
    loss_true = loss(P_xr_est, P_xr)
    cli_alert("Loss: {number(loss_true, 0.001)}")

    # from data
    cli_h2("From data")
    cli_ul()
    cli_li("{comma(N)} observations")
    cli_li("{number(N/n_x/n_s, 0.01)} average obs. per cell")
    cli_li("P_rs condition number: {number(kappa(P_rs), 0.1)}")
    cli_end()
    cli_text("\n")

    R = sample(n_r, N, replace=T, prob=P_r)
    P_sr = t(P_rs) %*% diag(1/P_r)
    S = samp_cond(N, R, P_sr)
    X = samp_cond(N, R, P_xr)
    if (impute) {
        N_est = 20
        R_ests = map(1:N_est, ~ samp_cond(N, S, P_rs %*% diag(1/colSums(P_rs))))
    }

    P_xs_est = prop.table(table(X, S))
    P_xr_act = prop.table(table(X, R), 2)

    R = sample(n_r, N, replace=T, prob=P_r)
    S = samp_cond(N, R, P_sr)
    P_rs_est = prop.table(table(R, S))

    P_xr_est1 = P_xs_est %*% MASS::ginv(P_rs)
    P_xr_est2 = P_xs_est %*% MASS::ginv(P_rs_est)

    P_xr_indep = matrix(rep(prop.table(table(X)), n_r), ncol=n_r)

    P_xr_est3 = P_xr_est1
    P_xr_est3[P_xr_est3 < 0] = 0
    P_xr_est3[P_xr_est3 > 1] = 1

    P_xr_est4 = 0.9*P_xr_est3 + 0.1*P_xr_indep

    if (impute) {
        P_xr_est5 = map(R_ests, function(R_est) prop.table(table(X, R_est), 2)) %>%
            reduce(~ .x + .y) %>%
            `/`(N_est)
    }

    loss_est1 = loss(P_xr_est1, P_xr)
    loss_est2 = loss(P_xr_est2, P_xr)
    loss_est3 = loss(P_xr_est3, P_xr)
    loss_est4 = loss(P_xr_est4, P_xr)
    if (impute) loss_est5 = loss(P_xr_est5, P_xr)
    loss_indep = loss(P_xr_indep, P_xr)
    cli_alert("Loss (independence): {number(loss_indep, 0.001)}")
    cli_alert("Rel. loss (true R,S): {number(loss_est1/loss_indep, 0.01)}")
    cli_alert("Rel. loss (corrected, true): {number(loss_est3/loss_indep, 0.01)}")
    cli_alert("Rel. loss (shrinkage, true): {number(loss_est4/loss_indep, 0.01)}")
    cli_alert("Rel. loss (est. R,S): {number(loss_est2/loss_indep, 0.01)}")
    if (impute) cli_alert("Rel. loss (imputation): {number(loss_est5/loss_indep, 0.01)}")
    invisible(NULL)
}

run_est(5, 40, 5e3, N=2.5e4, s_strength=0.8, loss=l_rmse)
