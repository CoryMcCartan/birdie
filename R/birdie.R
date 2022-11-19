#' @export
birdie <- function(r_probs, formula, data=NULL,
                   alpha=NULL, prefix="pr_", iter=20) {
    # if (missing(data)) cli_abort("{.arg data} must be provided.")
    # d_model = model.frame(formula, data=data, na.action=na.fail)
    # Y_vec = model.response(d_model)
    Y_vec = lme4::glFormula(formula, data=data)$fr[[1]]

    if (!is.matrix(r_probs)) {
        r_probs = as.matrix(dplyr::select(r_probs, starts_with(prefix)))
        colnames(r_probs) = substring(colnames(r_probs), nchar(prefix)+1L)
    }

    if (!check_vec(Y_vec)) cli_abort("{.arg X} must be a character or factor with no missing values.")
    if (nrow(data) != nrow(r_probs))
        cli_abort("{.arg data} and {.arg r_probs} must have the same number of rows.")

    Y_vec = as.factor(Y_vec)

    if (is.null(alpha)) {
        cli_inform("Using uniform prior for {.arg alpha} = Pr(X | R)")
        alpha = rep(1, nlevels(Y_vec))
    }
    if (length(alpha) != nlevels(Y_vec))
        cli_abort("{.arg alpha} prior must have the same number of elements as there are levels of X")



    out = em_nocov(as.integer(Y_vec), r_probs, alpha, iter=iter)
    # out = em_lmer(Y_vec, r_probs, formula, data, iter=iter)

    colnames(out$map) = colnames(r_probs)
    rownames(out$map) = levels(Y_vec)
    colnames(out$p_ryxs) = stringr::str_c(prefix, colnames(r_probs))
    out$p_ryxs = as_tibble(out$p_ryxs)
    out$map0 = NULL
    out$N = length(Y_vec)
    # out$vars = GZ_names
    out$x_lev = levels(Y_vec)
    out$r_lev = colnames(r_probs)
    class(out) = "birdie"

    out
}

em_lmer <- function(Y, p_rxs, form, data, iter=10, tol=0.001) {
    rlang::check_installed("lme4", "for this fitting method.")

    n_y = nlevels(Y)
    n_r = ncol(p_rxs)
    Y = as.integer(Y)
    # init with (penalized) weighted estimator
    p_ryxs = p_rxs
    est = dirichlet_map(Y, p_rxs, rep(2, n_y))
    p_ryxs = calc_bayes(Y, est, p_rxs)
    cat("E")

    idx_grp = lme4::glFormula(form, data=data)$fr[-1] |>
        dplyr::group_by_all() |>
        dplyr::group_indices()
    n_grp = max(idx_grp)
    grp_revlk = lapply(seq_len(n_grp), \(i) which(idx_grp == i))
    idx_extr = vapply(grp_revlk, \(x) x[1], 1L)

    ests = array(dim=c(n_grp, n_y, n_r))
    updating = matrix(TRUE, n_y, n_r)
    form_fit = update(form, cbind(succ, fail) ~ .)
    form_env = rlang::f_env(form_fit)
    break_next = FALSE
    for (i in seq_len(iter)) {
        last_ests = ests
        # M step
        for (r in seq_len(n_r)) {
            for (y in seq_len(n_y)) {
                if (!updating[y, r]) next

                rlang::env_bind(
                    form_env,
                    succ = tapply((Y == y) * p_ryxs[, r], idx_grp, sum),
                    fail = tapply((Y != y) * p_ryxs[, r], idx_grp, sum)
                )
                m = lme4::glmer(form_fit,
                                data = data[idx_extr, ],
                                family=binomial(), nAGQ=0L) |>
                    suppressWarnings() |>
                    suppressMessages()

                ests[, y, r] = fitted(m)
            }

            # normalize
            ests[, , r] = ests[, , r] / rowSums(ests[, , r])

            # check convergence
            if (i > 1) {
                for (y in seq_len(n_y)) {
                    rel_diff = mean(abs((ests[, y, r] - last_ests[, y, r]) / last_ests[, y, r]))
                    if (rel_diff <= tol) updating[y, r] = FALSE
                }
            }
        }
        cat("M")

        # E step
        for (j in seq_len(n_grp)) {
            idx = grp_revlk[[j]]
            p_ryxs[idx, ] = calc_bayes(Y[idx], ests[j, , ],
                                       p_rxs[idx, , drop=FALSE])
        }
        cat("E")

        if (break_next) break
        if (all(!updating)) { # converged -- do one more full update
            updating = matrix(TRUE, n_y, n_r)
            break_next = TRUE
        }
    }
    cat("\n")

    # final global mean
    est = dirichlet_map(Y, p_ryxs, rep(1, n_y))

    list(map = est,
         ests = ests,
         p_ryxs = p_ryxs)
}


#' @export
print.birdie <- function(x, ...) {
    cli::cli_text("A {.pkg BIRDiE} model fit with
                  {format(x$N, big.mark=',')} observations")
    cat("\n")

    cli::cli_text("Estimated distribution:")
    cat("\n")
    m = round(x$map, 3)
    colnames(m) = x$r_lev
    rownames(m) = x$x_lev
    print(m)
}

