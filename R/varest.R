
# bootstrap `em_cat_dir()`
boot_cat_dir <- function(mle, R=10, Y, X, weights, p_rxs, prior, n_x, ctrl) {
    N = length(Y)
    W_tot = sum(weights)
    n_r = ncol(prior$alpha)
    n_y = nrow(prior$alpha)

    out = matrix(nrow=length(prior$alpha), ncol=R)

    ctrl$abstol = max(0.0005, ctrl$abstol)
    ctrl$reltol = max(0.005, ctrl$reltol)
    ctrl$max_iter = 50

    ones = rep_along(Y, 1)
    ones_mat = matrix(1, nrow=n_y, ncol=n_r)
    mk_wt = weight_maker(N, R, weights)

    cli::cli_progress_bar("Bootstrapping", total=R)
    for (i in seq_len(R)) {
        wt = mk_wt()

        res = ctrl$accel(mle, function(curr) {
            # reproject if acceleration has brought us out of bounds
            curr[curr < 0] = 0 + 1e3*.Machine$double.eps
            curr[curr > 1] = 1 - 1e3*.Machine$double.eps

            .Call(`_birdie_em_dirichlet_wt`, curr, Y, X, wt, p_rxs, prior$alpha, n_x)
        }, ctrl, n_x=n_x)

        out[, i] = em_dirichlet_wt(res$ests, Y, ones, wt, p_rxs, ones_mat, 1)

        cli::cli_progress_update()
    }
    cli::cli_progress_done()

    out
}

boot_lm <- function(mle, R=10, Y, X, weights, p_rxs, prior, ctrl) {
    N = length(Y)
    p = ncol(X)
    ign = -seq_len(p)
    n_r = ncol(p_rxs)
    W_tot = sum(weights)

    out = matrix(nrow=n_r, ncol=R)

    ctrl$abstol = max(0.0005, ctrl$abstol)
    ctrl$reltol = max(0.005, ctrl$reltol)
    ctrl$max_iter = 50

    mk_wt = weight_maker(N, R, weights)

    cli::cli_progress_bar("Bootstrapping", total=R)
    for (i in seq_len(R)) {
        wt = mk_wt()

        res = ctrl$accel(mle, function(curr) {
            coefs = matrix(curr[-1], nrow=p, ncol=n_r)
            p_ryxs = lm_estep(X, Y, coefs, curr[1], p_rxs, p)
            lm_mstep(X, Y, p_ryxs, wt, TRUE, p, prior, ctrl)
        }, ctrl, n_x=p)

        ests = matrix(res$ests[-1], nrow=p, ncol=n_r)
        p_ryxs = lm_estep(X, Y, ests, res$ests[1], p_rxs, p)
        out[, i] = colSums(p_ryxs * (X[ign, ] %*% ests)) / colSums(p_ryxs)

        cli::cli_progress_update()
    }
    cli::cli_progress_done()

    out
}


weight_maker <- function(N, R, weights) {
    if (N > 1000 && R > 100) {
        function() tabulate(sample.int(N, sum(weights), replace=TRUE), N) / N
    } else { # more computationally intensive but smoother
        function() as.numeric(rdirichlet(1, weights))
    }
}
