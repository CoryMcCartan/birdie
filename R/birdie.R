#' @export
birdie <- function(r_probs, formula, data=NULL,
                   prior=NULL, prefix="pr_", max_iter=30) {
    # if (missing(data)) cli_abort("{.arg data} must be provided.")

    # figure out type of model and extract response vector
    re_terms = lme4::findbars(formula)
    if (is.null(re_terms)) {
        d_model = model.frame(formula, data=data, na.action=na.fail)
        Y_vec = model.response(d_model)
        method = "nocov"
    } else {
        d_model = lme4::lFormula(formula, data=data)$fr
        Y_vec = d_model[[1]]
        method = if (length(re_terms) == 1) "re1" else "glmm"
    }

    # set up race probability matrix
    if (!is.matrix(r_probs)) {
        r_probs = as.matrix(select(r_probs, starts_with(prefix)))
        colnames(r_probs) = substring(colnames(r_probs), nchar(prefix)+1L)
    }
    # check predictors
    if (inherits(r_probs, "bisg")) {
        if (attr(r_probs, "S_name") %in% colnames(d_model)) {
            cli_warn("Last name vector {.arg {attr(r_probs, 'S_name')}}
                     should not be used in BIRDiE model")
        }
        if (length(setdiff(colnames(d_model), attr(r_probs, "GX_names"))) > 0) {
            cli_warn(c("Found variables used in BIRDiE model which were not
                     used to create BISG probabilities.",
                     "x"="Statistically valid inference is not guaranteed."))
        }
    }

    # check types
    if (!check_vec(Y_vec)) cli_abort("{.arg X} must be a character or factor with no missing values.")
    if (nrow(data) != nrow(r_probs))
        cli_abort("{.arg data} and {.arg r_probs} must have the same number of rows.")

    Y_vec = as.factor(Y_vec)

    if (is.null(prior)) {
        if (method == "nocov")
            cli_inform("Using uniform prior for {.arg alpha} = Pr(X | R)")
        prior = rep(1, nlevels(Y_vec))
    }
    if (length(prior) != nlevels(Y_vec))
        cli_abort("{.arg alpha} prior must have the same number of elements as there are levels of X")

    # run inference
    if (method == "nocov") {
        out = em_nocov(as.integer(Y_vec), r_probs, prior, iter=max_iter)
    } else if (method == "re1") {
        out = em_re1(Y_vec, r_probs, formula, data, iter=max_iter)
    } else if (method == "lmer") {
        out = em_glmm(Y_vec, r_probs, formula, data, iter=max_iter)
    }


    # output
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

em_re1 <- function(Y, p_rxs, form, d_model, prior=rep(1, ncol(p_rxs)),
                   iter=10, tol=0.001) {
    n_y = nlevels(Y)
    n_r = ncol(p_rxs)
    Y = as.integer(Y)
    # init with (penalized) weighted estimator
    p_ryxs = p_rxs
    est = dirichlet_map(Y, p_rxs, rep(2, n_y))
    p_ryxs = calc_bayes(Y, est, p_rxs)

    d_model =
    idx_grp = vctrs::vec_duplicate_id(lme4::lFormula(form, data=data)$fr[-1]) |>
        as.factor() |>
        as.integer()
    n_grp = max(idx_grp)
    grp_revlk = lapply(seq_len(n_grp), \(i) which(idx_grp == i))
    idx_extr = vapply(grp_revlk, \(x) x[1], 1L)

    ests = array(dim=c(n_grp, n_y, n_r))

    form_fit = update.formula(form, cbind(succ, fail) ~ .)
    form_env = rlang::f_env(form_fit)

    break_next = FALSE
    for (i in seq_len(iter)) {
        last_ests = ests
        # M step
        for (r in seq_len(n_r)) {

            rlang::env_bind(
                form_env,
                succ = sum_grp((Y == y) * p_ryxs[, r], idx_grp, n_grp),
                fail = sum_grp((Y != y) * p_ryxs[, r], idx_grp, n_grp)
            )
            m = lme4::glmer(form_fit,
                            data = data[idx_extr, ],
                            family=binomial(), nAGQ=0L) |>
                suppressWarnings() |>
                suppressMessages()

            ests[, y, r] = fitted(m)

        }

        # E step
        for (j in seq_len(n_grp)) {
            idx = grp_revlk[[j]]
            p_ryxs[idx, ] = calc_bayes(Y[idx], ests[j, , ],
                                       p_rxs[idx, , drop=FALSE])
        }

        if (i > 1) {
            do_break = FALSE
            for (y in seq_len(n_y)) {
                rel_diff = max(abs((ests[, y, r] - last_ests[, y, r]) /
                                       last_ests[, y, r]))
                if (rel_diff <= tol) do_break = TRUE
            }
            if (do_break) break
        }
    }

    # final global mean
    est = dirichlet_map(Y, p_ryxs, rep(1, n_y))

    list(map = est,
         ests = ests,
         p_ryxs = p_ryxs)
}

em_glmm <- function(Y, p_rxs, form, data, iter=10, tol=0.001) {
    n_y = nlevels(Y)
    n_r = ncol(p_rxs)
    Y = as.integer(Y)
    # init with (penalized) weighted estimator
    p_ryxs = p_rxs
    est = dirichlet_map(Y, p_rxs, rep(2, n_y))
    p_ryxs = calc_bayes(Y, est, p_rxs)

    idx_grp = vctrs::vec_duplicate_id(lme4::lFormula(form, data=data)$fr[-1]) |>
        as.factor() |>
        as.integer()
    n_grp = max(idx_grp)
    grp_revlk = lapply(seq_len(n_grp), \(i) which(idx_grp == i))
    idx_extr = vapply(grp_revlk, \(x) x[1], 1L)

    ests = array(dim=c(n_grp, n_y, n_r))
    updating = matrix(TRUE, n_y, n_r)

    form_fit = update.formula(form, cbind(succ, fail) ~ .)
    form_env = rlang::f_env(form_fit)

    break_next = FALSE
    cli::cli_progress_bar("Fitting BIRDiE with EM", total=iter*n_y*n_r)
    for (i in seq_len(iter)) {
        last_ests = ests
        # M step
        for (r in seq_len(n_r)) {
            for (y in seq_len(n_y)) {
                cli::cli_progress_update()
                if (!updating[y, r]) next

                rlang::env_bind(
                    form_env,
                    succ = sum_grp((Y == y) * p_ryxs[, r], idx_grp, n_grp),
                    fail = sum_grp((Y != y) * p_ryxs[, r], idx_grp, n_grp)
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

        # E step
        for (j in seq_len(n_grp)) {
            idx = grp_revlk[[j]]
            p_ryxs[idx, ] = calc_bayes(Y[idx], ests[j, , ],
                                       p_rxs[idx, , drop=FALSE])
        }

        if (break_next) break
        if (all(!updating)) { # converged -- do one more full update
            updating = matrix(TRUE, n_y, n_r)
            break_next = TRUE
        }
    }
    cli::cli_progress_done()

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

