
#' Estimate X versus Race Given Name, Location, Outcome, and Covariates
#'
#' Bayesian model for X|R given name and covariates
#'
#' @param r_probs a matrix data frame of race probabilities
#' @param X the column containing the outcome
#' @param G the column containing locations
#' @param Z the column(s), if any, containing other covariates. Use `c()` to provide multiple columns
#' @param data the data
#' @param prefix how to select the race probability columns from `r_probs`, if
#'   if is a data frame
#' @param config a list containing extra arguments to the Pyro model
#' @param silent whether to suppress output
#' @param reload_py whether to reload the python code (development only)
#'
#' @return TBD
#' @export
model_race = function(r_probs, X, G, Z=NULL, data=NULL, sgz=NULL, prefix="pr_",
                      config=list(), silent=FALSE, reload_py=FALSE) {
    if (missing(data)) cli_abort("{.arg data} must be provided.")
    X_vec = eval_tidy(enquo(X), data)
    G_vec = eval_tidy(enquo(G), data)
    Z_df = data[, eval_select(enquo(Z), data)]
    if (!is.matrix(r_probs)) {
        r_probs = as.matrix(dplyr::select(r_probs, starts_with(prefix)))
        colnames(r_probs) = substring(colnames(r_probs), nchar(prefix)+1L)
    }

    if (!check_vec(X_vec)) cli_abort("{.arg X} must be a character or factor with no missing values.")
    if (nrow(data) != nrow(r_probs))
        cli_abort("{.arg data} and {.arg r_probs} must have the same number of rows.")
    if (!is.character(G_vec) && !is.factor(G_vec))
        cli_abort("{.arg G} must be a character or factor vector.")
    if (!all(vapply(Z_df, class, character(1)) %in% c("character", "factor"))) {
        cli_abort("{.arg Z} must contain only character or factor columns.")
    }
    if (any(is.na(Z_df))) cli_abort("Missing values found in {.arg Z}")

    X_vec = as.factor(X_vec)
    G_vec = coalesce(G_vec, "<none>")
    GZ = cbind(as.factor(G_vec), Z_df)
    GZ_levels = vapply(GZ, nlevels, integer(1))
    # set up matrix
    GZ_var = inverse.rle(list(lengths=GZ_levels, values=1:ncol(GZ)))
    suppressWarnings({
        GZ_mat = do.call(cbind, lapply(GZ, function(x) {
            out = matrix(0, nrow=length(x), ncol=nlevels(x))
            out[cbind(seq_along(x), as.integer(x))] = 1
            out
        }))
    })

    if (isTRUE(reload_py) || !"fit_additive" %in% names(py)) {
        reticulate::py_run_string("if 'py.utils' in sys.modules.keys(): del sys.modules['py.utils']")
        reticulate::py_run_string("if 'py.fit' in sys.modules.keys(): del sys.modules['py.fit']")
        reticulate::py_run_string("from tqdm import tqdm; tqdm._instances.clear()")
        reticulate::py_run_file(system.file("py/pyro.py", package="raceproxy"))
    }

    defaults = list(max_iter = 5000,
                    subsamp = 2048,
                    epoch = 50,
                    draws = 800,
                    lr = 0.2,
                    n_mi = 5,
                    it_avgs = 300,
                    tol_rhat = 1.2)
    for (i in names(defaults)) {
        if (is.null(config[[i]]))
            config[[i]] = defaults[[i]]
    }

    prior = list(x = 5.00, xr = 0.75, beta = 1.00)

    sgz$pr = NULL
    #sgz$n_s = nlevels(sgz$S)
    sgz$n_gz = nlevels(sgz$GZ)
    #sgz$S = as.integer(sgz$S)
    sgz$GZ = as.integer(sgz$GZ)

    tictoc::tic()
    out = py$fit_additive(as.integer(X_vec), GZ_mat, as.integer(GZ_var), r_probs,
                          nlevels(X_vec), max(GZ_var), sgz, prior,
                          it=as.integer(config$max_iter),
                          epoch=as.integer(config$epoch),
                          subsamp=min(length(X_vec), as.integer(config$subsamp)),
                          n_draws=as.integer(config$draws),
                          it_avgs=as.integer(config$it_avgs),
                          n_mi=as.integer(config$n_mi),
                          lr=config$lr, tol_rhat=config$tol_rhat,
                          silent=silent)
    tictoc::toc(quiet=silent)

    out
}
