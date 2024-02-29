# Minimal functions to fit Stan model
# Assumes everything is correctly formatted -- all safeguards turned off
# Good general guide here <https://cran.r-project.org/web/packages/StanHeaders/vignettes/stanmath.html>
# Also see rstan optimizing code <https://github.com/stan-dev/rstan/blob/develop/rstan/rstan/R/stanmodel-class.R#L346>
# This file is called stanmodels.R so that rstantools::rstan_config() won't create one

get_stanmodel <- function(module, data, seed=5118L) {
    methods::new(module, data, seed, function() stop("cxxfn"))
}

get_skeleton <- function(mod) {
    par_nm = mod$param_names()
    p_dims = mod$param_dims()[par_nm != "lp__"]
    lapply(p_dims, function(d) array(0, dim=d))
}

optim_model_stan <- function(mod, init, skeleton,
                             algorithm = c("LBFGS"), ..., seed=5118L) {
    if (is.numeric(init))
        init = as.character(init)

    args = list(init = init,
                seed = as.integer(seed),
                method = "optim",
                refresh = 0L,
                algorithm = match.arg(algorithm),
                ...)

    res = mod$call_sampler(args)

    theta = relist(res$par, skeleton)

    list(
        par = mod$unconstrain_pars(theta),
        converged = (attr(res, "return_code") == 0)
    )
}

optim_model <- function(mod, init, tol_rel_obj, ...) {
    if (is.numeric(init))
        init = as.character(init)

    res = optim(init,
                fn = \(x) -mod$log_prob(x, TRUE, FALSE),
                gr = \(x) -mod$grad_log_prob(x, FALSE),
                method = "L-BFGS-B",
                control = list(factr=0.1/tol_rel_obj, maxit=500))

    list(
        par = res$par,
        converged = (res$convergence == 0)
    )
}

constrain_pars <- function(mod, skeleton, x) {
    relist(mod$constrain_pars(x), skeleton)
}
