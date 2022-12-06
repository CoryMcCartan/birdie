# EM Acceleration

accel_none <- function(init, em_step, ctrl, n_x=1, ...) {
    ests = init

    converged = FALSE
    incr_factor = n_x^(1/2)
    for (i in seq_len(ctrl$max_iter)) {
        last = ests
        ests = em_step(ests)

        if (check_convergence(ests, last,
                              ctrl$abstol*incr_factor, ctrl$reltol*incr_factor)) {
            converged = TRUE
            break
        }
    }

    absdiff = abs(ests - last)
    list(ests = ests,
         absdiff = max(absdiff),
         reldiff = max(absdiff / ests),
         iters = i,
         converge = converged)
}


accel_anderson <- function(init, em_step, ctrl, n_x=1, ...) {
    if (ctrl$order <= 0) ctrl$order = min(12, length(init)/2)

    # iterations matrix
    X = cbind(em_step(init), init)
    # residuals matrix
    f = cbind(em_step(X[, 1]), X[, 1]) - X
    # residual increments matrix (first element telescopes)
    dF = f[, 1, drop=FALSE] - f[, 2, drop=FALSE]
    # increments matrix
    dX = X[, 1, drop=FALSE] - X[, 2, drop=FALSE]

    converged = FALSE
    incr_factor = n_x^(1/2)
    # dampen = 1e-12
    for (i in seq_len(ctrl$max_iter)[-1]) { # i = 2 ... max
        # predict residuals from residual increments
        k = ncol(dF)
        # fit penalized least squares
        # gamma = .lm.fit(rbind(dF, diag(k)*dampen), c(f[, 1], rep_len(0, k)))$coefficients
        gamma = .lm.fit(dF, f[, 1])$coefficients

        # do Anderson step
        new_X = X[, 1] + f[, 1] -  (dX + dF) %*% gamma
        new_X[new_X < 0] = 0
        new_X[new_X > 1] = 1
        X = cbind(new_X, X)
        f = cbind(em_step(X[, 1]) - X[, 1], f)

        # update increments
        dF = cbind(f[, 1] - f[, 2], dF)
        dX = cbind(X[, 1] - X[, 2], dX)

        m_k = min(ctrl$order, i)
        if (ncol(dX) > ctrl$order) {
            if (ctrl$restart) { # restart
                X = X[, 1:2]
                f = f[, 1:2]
                dX = dX[, 1, drop=FALSE]
                dF = dF[, 1, drop=FALSE]
            } else { # truncate to order
                cols = seq_len(ctrl$order + 1)
                X = X[, cols]
                f = f[, cols]
                cols = seq_len(ctrl$order)
                dX = dX[, cols]
                dF = dF[, cols]
            }
        }

        if (check_convergence(X[, 1], X[, 2],
                              ctrl$abstol*incr_factor, ctrl$reltol*incr_factor)) {
            converged = TRUE
            break
        }
    }

    ests = X[, 1] + f[, 1]
    absdiff = ests - X[, 2] - f[, 2]
    list(ests = X[, 1] + f[, 1], # use extra em_step we already computed
         absdiff = max(absdiff),
         reldiff = max(absdiff / ests),
         iters = i,
         converge = converged)
}

accel_squarem <- function(init, em_step, ctrl, n_x=1, ...) {
    # rlang::check_installed("SQUAREM", "For SQUAREM acceleration.")
    squarem_method = if (ctrl$order == 1) 3 else "rre"
    incr_factor = n_x^(1/2)

    res = SQUAREM::squarem(
        init,
        em_step,
        control=list(K=ctrl$order, method=squarem_method, minimize=TRUE,
                     square=TRUE, step.min0=0.5, step.max0=1, mstep=4, objfn.inc=1,
                     kr=1, tol=ctrl$abstol*incr_factor, maxiter=ctrl$max_iter,
                     trace=FALSE, intermed=FALSE)
    )

    list(ests = res$par,
         iters = res$fpeval,
         converge = res$convergence)
}


accel_daarem <- function(init, em_step, ctrl, n_x=1, ...) {
    rlang::check_installed("daarem", "For DAAREM acceleration.")
    if (ctrl$order <= 0) ctrl$order = min(10, length(init)/2)
    incr_factor = n_x^(1/2)

    res = suppressWarnings(daarem::daarem(
        init,
        em_step,
        control=list(maxiter=ctrl$max_iter, order=ctrl$order,
                     tol=ctrl$abstol*incr_factor, mon.tol=0.01, cycl.mon.tol=0.0,
                     alpha=1.2, kappa=25, resid.tol=0.95, convtype="param")
    ))

    list(ests = res$par,
         iters = res$fpeval,
         converge = res$convergence)
}
