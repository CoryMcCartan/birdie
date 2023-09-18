
# bootstrap `em_cat_dir()`
boot_cat_dir <- function(mle, R=10, Y, X, p_rxs, prior, n_x, ctrl) {
    N = length(Y)
    n_r = ncol(prior$alpha)
    n_y = nrow(prior$alpha)

    out = matrix(nrow=length(prior$alpha), ncol=R)

    ctrl$abstol = 0.0005 #min(ctrl$abstol * 1000, 0.001)
    ctrl$reltol = 0.005 #min(ctrl$reltol * 1000, 0.001)
    ctrl$max_iter = 50

    ones = rep_along(Y, 1)
    ones_mat = matrix(1, nrow=n_y, ncol=n_r)

    if (N > 1000 && R > 100) {
        mk_wt = function() tabulate(sample.int(N, N, replace=TRUE), N) / N
    } else { # more computationally intensive but smoother
        mk_wt = function() as.numeric(rdirichlet(1, ones))
    }

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
