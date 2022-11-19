
#' Estimate X versus Race Given Name, Location, Outcome, and Covariates
#'
#' Bayesian model for X|R given name and covariates
#'
#' @param r_probs a matrix data frame of race probabilities
#' @param X the column containing the outcome
#' @param G the column containing locations
#' @param Z the column(s), if any, containing other covariates. Use `c()` to provide multiple columns.
#' @param condition the columns in `G` or `Z` to condition on in creating predictions.
#'   For example, if `Z` contains `gender`, setting `condition=gender` would
#'   produce X|R estimates for each gender.
#' @param data the data
#' @param prefix how to select the race probability columns from `r_probs`, if
#'   if is a data frame
#' @param config a list containing extra arguments to the Pyro model
#' @param silent whether to suppress output
#' @param reload_py whether to reload the python code (development only)
#'
#' @return A list containing the model output. Element `p_xr` contains the
#'   approximate posterior draws of the global X|R table.
#' @noRd
model_race = function(r_probs, X, G, Z=NULL, condition=NULL,
                      data=NULL, prefix="pr_", method=c("svi", "mle"),
                      config=list(), silent=FALSE, reload_py=FALSE) {
    if (missing(data)) cli_abort("{.arg data} must be provided.")
    X_vec = eval_tidy(enquo(X), data)
    G_vec = eval_tidy(enquo(G), data)
    Z_df = data[, tidyselect::eval_select(enquo(Z), data)]
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
            out = matrix(0L, nrow=length(x), ncol=nlevels(x))
            out[cbind(seq_along(x), as.integer(x))] = 1L
            out
        }))
    })


    # prepare prediction matrices
    GZ_names = c(
        names(tidyselect::eval_select(enquo(G), data)),
        names(tidyselect::eval_select(enquo(Z), data))
    )
    cond_name = names(tidyselect::eval_select(enquo(condition), data))

    preds = list(global = colMeans(GZ_mat))
    if (length(cond_name) == 1) {
        col_idx = match(cond_name, GZ_names)
        cols = which(GZ_var == col_idx)
        col_lvls = levels(GZ[[col_idx]])
        for (i in seq_along(cols)) {
            lab = paste0(cond_name, ": ", col_lvls[i])
            preds[[lab]] = colMeans(GZ_mat[GZ_mat[, cols[i]] == 1, , drop=FALSE])
        }
    }

    method = match.arg(method)

    defaults = list(max_iter = 5000,
                    subsamp = 2048,
                    epoch = 50,
                    draws = if_else(method == "svi", 1000, 2),
                    lr = if_else(method == "svi", 0.25, 0.3),
                    n_err = 100,
                    prior = list(x = 5.00, xr = 0.75, beta = 1.00),
                    it_avgs = if_else(method == "svi", 400, 600),
                    tol_rhat = if_else(method == "svi", 1.2, 1.1))
    for (i in names(defaults)) {
        if (is.null(config[[i]]))
            config[[i]] = defaults[[i]]
    }


    ts1 = proc.time()
    out = list()

    class(out) = "fit_birdie"
    out$N = length(X_vec)
    out$vars = GZ_names
    out$x_lev = levels(X_vec)
    out$r_lev = colnames(r_probs)
    out
}

#' @export
print.fit_birdie = function(x, ...) {
    cli::cli_text("A {.pkg BIRDiE} model fit with
                  {format(x$N, big.mark=',')} observations and
                  {format(dim(x$draws$global)[1], big.mark=',')} draws")
    # cli::cli_text("{dim(fit$draws$global)[2]} outcome and
                  # {dim(fit$draws$global)[3]} race categories")
    cat("\n")

    cli::cli_text("Predictor variable importance (std. devs.):")
    if (is.matrix(x$beta_scale)) {
        var_sd = round(colMeans(x$beta_scale), 3)
    } else {
        var_sd = mean(x$beta_scale)
    }
    names(var_sd) = x$vars
    print(var_sd)
    cat("\n")

    cli::cli_text("Predictions for: {c('everyone', names(x$draws)[-1])}")
    cat("\n")

    cli::cli_text("Estimates for everyone:")
    m = round(colMeans(x$draws$global), 3)
    colnames(m) = x$r_lev
    rownames(m) = x$x_lev
    print(m)
}

