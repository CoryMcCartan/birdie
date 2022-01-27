est_nonparam = function(X, GZ, pr_base, alpha, method, iter, verbose) {
    stan_data = list(
        N = length(X),
        n_x = nlevels(X),
        n_r = length(p_r),
        n_gz = nlevels(GZ),

        X = as.integer(X),
        GZ = as.integer(GZ),
        pr_base = pr_base,
        p_gz = prop.table(table(GZ)),

        n_prior_obs = alpha[1]
    )

    tictoc::tic()
    if (method == "opt") {
        out = rstan::optimizing(stanmodels$simple_nonparam, stan_data,
                                verbose=verbose,
                                #draws=iter, importance_resampling=TRUE,
                                tol_grad=1e-6, tol_param=1e-5, init_alpha=1)
    } else if (method == "vb") {
        out = rstan::vb(stanmodels$simple_nonparam, stan_data,
                        init=0,
                        pars="p_xr", algorithm="meanfield",
                        eta=1, adapt_engaged=FALSE,
                        grad_samples=2, eval_elbo=50, elbo_samples=75,
                        iter=700, importance_resampling=FALSE)
    } else {
        out = rstan::sampling(stanmodels$simple_nonparam, stan_data, pars="p_xr",
                              chains=1, iter=300+iter, warmup=300)
    }
    tictoc::toc()
    out
}

est_additive = function(X, GZ, GZ_var, pr_base, method, iter, verbose) {
    stan_data = list(
        N = length(X),
        n_x = nlevels(X),
        n_r = length(p_r),
        n_gz = max(GZ_var),
        n_gz_col = ncol(GZ),

        X = as.integer(X),
        GZ = GZ,
        n_nonzero = as.integer(sum(GZ)),
        GZ_var = GZ_var,
        pr_base = pr_base,

        prior_lp_x_scale = 5,
        prior_lp_xr_scale = 0.75,
        prior_beta_scale = 1.0
    )

    tictoc::tic()
    if (method == "opt") {
        out = rstan::optimizing(stanmodels$simple_additive, stan_data,
                                verbose=verbose,
                                #draws=iter, importance_resampling=TRUE,
                                tol_grad=1e-6, tol_param=1e-5, init_alpha=1)
    } else if (method == "vb") {
        out = rstan::vb(stanmodels$simple_additive, stan_data,
                        init=0,
                        pars=c("p_xr", "beta_scale"), algorithm="meanfield",
                        eta=1, adapt_engaged=FALSE,
                        grad_samples=1, eval_elbo=50, elbo_samples=75,
                        tol_rel_obj=0.002,
                        iter=700, importance_resampling=FALSE)
    } else {
        out = rstan::sampling(stanmodels$simple_additive, stan_data, pars="p_xr",
                              #init=list(list(p_xr=xr$bisg)),
                              chains=1, iter=300+iter, warmup=300)
    }
    tictoc::toc()
    out
}
