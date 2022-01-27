est_nonparam_pyro = function(X, GZ, pr_base, alpha,
                             iter=2000, subsamp=1000, smooth=200, lr=0.01, tol=0.005,
                             reload=F) {
    if (isTRUE(reload)) {
        reticulate::py_run_string("if 'py.utils' in sys.modules.keys(): del sys.modules['py.utils']")
        reticulate::py_run_file(system.file("py/pyro.py", package="raceproxy"))
    }

    p_gz = as.array(prop.table(table(GZ)))

    tictoc::tic()
    out = py$fit_nonparam(as.integer(X), as.integer(GZ), pr_base,
                          p_gz, nlevels(X), alpha,
                          as.integer(iter), as.integer(subsamp),
                          lr, as.integer(smooth), tol)
    tictoc::toc()
    out
}

est_additive_pyro = function(X, GZ, GZ_var, pr_base,
                             iter=2000, subsamp=1000, lr=0.01, tol_rhat=1.3,
                             reload=F) {
    if (isTRUE(reload)) {
        reticulate::py_run_string("if 'py.utils' in sys.modules.keys(): del sys.modules['py.utils']")
        reticulate::py_run_file(system.file("py/pyro.py", package="raceproxy"))
    }

    prior = list(x = 5, xr = 0.75, beta = 1.0)

    tictoc::tic()
    out = py$fit_additive(as.integer(X), GZ, as.integer(GZ_var), pr_base,
                          nlevels(X), max(GZ_var), prior,
                          it=as.integer(iter), subsamp=as.integer(subsamp),
                          lr=lr, tol_rhat=tol_rhat)
    tictoc::toc()
    out
}
