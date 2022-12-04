# Minimal functions to fit Stan model
# Assumes everything is correctly formatted -- all safeguards turned off
# Good general guide here <https://cran.r-project.org/web/packages/StanHeaders/vignettes/stanmath.html>
# Also see rstan optimizing code <https://github.com/stan-dev/rstan/blob/develop/rstan/rstan/R/stanmodel-class.R#L346>

get_stanmodel <- function(module, data, seed=5118L) {
    new(module, data, seed, function() stop("cxxfn"))
}

get_skeleton <- function(mod) {
    par_nm = mod$param_names()
    p_dims = mod$param_dims()[par_nm != "lp__"]
    lapply(p_dims, function(d) array(0, dim=d))
}

optim_model <- function(mod, init, skeleton,
                        algorithm = c("LBFGS", "BFGS", "Newton"), ..., seed=5118L) {
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

constrain_pars <- function(mod, skeleton, x) {
    relist(mod$constrain_pars(x), skeleton)
}
