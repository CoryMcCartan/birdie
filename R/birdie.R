#' @export
birdie <- function(r_probs, formula, data=NULL,
                   prior=NULL, prefix="pr_", max_iter=30) {
    # if (missing(data)) cli_abort("{.arg data} must be provided.")

    # figure out type of model and extract response vector
    re_terms = lme4::findbars(formula)
    if (is.null(re_terms)) {
        d_model = model.frame(formula, data=data, na.action=na.fail)
        Y_vec = model.response(d_model)
        if (length(attr(terms(formula), "term.labels")) == 0) { # just Y ~ 1
            method = "pool"
        } else { # saturated
            method = "sat"

            tt = terms(formula)
            int_ord = attr(tt, "order")
            int_ord = table(factor(int_ord, levels=seq_len(max(int_ord))))
            full_int = all(int_ord == choose(length(int_ord), seq_along(int_ord))) # pascal's triangle
            if (!full_int) {
                x_vars = attr(tt, "term.labels")[attr(tt, "order") == 1]
                cli_warn(c("Saturated (no-pooling) model being fit without
                         full interaction structure. No estimation guarantees.",
                         "i"="For full interaction structure, use
                         `{colnames(d_model)[1]} ~ {paste0(x_vars, collapse=' * ')}`."))
            }
        }
    } else {
        d_model = lme4::lFormula(formula, data=data)$fr
        Y_vec = d_model[[1]]
        method = if (length(re_terms) == 1) "re1" else "glmm"
    }

    # check predictors
    if (inherits(r_probs, "bisg")) {
        if (attr(r_probs, "S_name") %in% colnames(d_model)) {
            cli_warn("Last name vector {.arg {attr(r_probs, 'S_name')}}
                     should not be used in BIRDiE model")
        }
        if (method != "pool" &&
            length(setdiff( attr(r_probs, "GX_names"), colnames(d_model) )) > 0) {
            cli_warn(c("Some variables used to create BISG probabilities are
                       not used in BIRDiE model.",
                     "x"="Statistically valid inference is not guaranteed."))
        }
    }
    # set up race probability matrix
    if (!is.matrix(r_probs)) {
        r_probs = as.matrix(select(r_probs, starts_with(prefix)))
        colnames(r_probs) = substring(colnames(r_probs), nchar(prefix)+1L)
    }

    # check types
    if (!check_vec(Y_vec)) cli_abort("{.arg X} must be a character or factor with no missing values.")
    if (nrow(data) != nrow(r_probs))
        cli_abort("{.arg data} and {.arg r_probs} must have the same number of rows.")

    Y_vec = as.factor(Y_vec)

    if (is.null(prior)) {
        prior = rep_len(1, nlevels(Y_vec))
        if (method == "pool") {
            cli_inform("Using uniform prior for {.arg alpha} = Pr(X | R)")
        } else if (method == "sat") {
            cli_inform("Using c(2, 2, ... 2) prior for {.arg alpha} = Pr(X | R)")
            prior = rep_len(2, nlevels(Y_vec))
        }
    }
    if (length(prior) != nlevels(Y_vec))
        cli_abort("{.arg alpha} prior must have the same number of elements as there are levels of X")

    # run inference
    if (method == "pool") {
        out = em_pool(as.integer(Y_vec), r_probs, prior, iter=max_iter)
    } else if (method == "sat") {
        X_vec = to_unique_ids(d_model[-1])
        n_x = max(X_vec)
        out = em_sat(as.integer(Y_vec), X_vec, r_probs, prior, n_x, iter=max_iter)
    } else if (method == "re1") {
        # out = em_re1(Y_vec, r_probs, formula, data, iter=max_iter)
        out = em_glmm(Y_vec, r_probs, formula, data, iter=max_iter)
    } else if (method == "lmer") {
        out = em_glmm(Y_vec, r_probs, formula, data, iter=max_iter)
    }

    if (isTRUE(out$max_iters)) {
        cli_warn(c("Maximum number of EM iterations hit.",
                   "i"="Estimates may not have converged.",
                   ">"="Consider increasing {.arg max_iter}."),
                 call=parent.frame())
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
    out$method = method
    class(out) = "birdie"

    out
}

em_re1 <- function(Y, p_rxs, form, d_model, prior=rep(1, ncol(p_rxs)),
                   iter=10, tol=0.001) {

}

em_glmm <- function(Y, p_rxs, form, data, iter=10, tol=0.001, ref_grp=1L) {
    n_y = nlevels(Y)
    n_r = ncol(p_rxs)
    Y = as.integer(Y)
    ones = rep_along(Y, 1)
    # init with (penalized) weighted estimator
    p_ryxs = p_rxs
    est = dirichlet_map(Y, ones, p_rxs, rep(2, n_y), 1)[[1]]
    p_ryxs = calc_bayes(Y, ones, list(est), p_rxs, 1)

    idx_grp = to_unique_ids(lme4::lFormula(form, data=data)$fr[-1])
    n_grp = max(idx_grp)
    grp_revlk = lapply(seq_len(n_grp), \(i) which(idx_grp == i))
    idx_extr = vapply(grp_revlk, function(x) x[1], 1L)
    idx_extr_y = lapply(seq_len(n_y), function(y) {
        intersect(idx_extr, which(Y %in% c(y, ref_grp)))
    })

    ests = array(dim=c(n_grp, n_y, n_r))
    updating = matrix(TRUE, n_y, n_r)

    form_fit = update.formula(form, cbind(succ, fail) ~ .)
    form_env = rlang::f_env(form_fit)

    break_next = FALSE
    cli::cli_progress_bar("Fitting BIRDiE with EM", total=iter*(n_y-1)*n_r)
    for (i in seq_len(iter)) {
        last_ests = ests
        # M step
        for (r in seq_len(n_r)) {
            cts = sum_grp(Y, idx_grp, p_ryxs[, r], rep_len(0, n_y), n_y, n_grp)

            for (y in seq_len(n_y)[-ref_grp]) {
                cli::cli_progress_update()
                if (!updating[y, r]) next

                rlang::env_bind(
                    form_env,
                    succ = cts[, y],
                    fail = cts[, ref_grp],
                )
                m = lme4::glmer(form_fit,
                                data = data[idx_extr, ],
                                family=binomial(), nAGQ=0L) |>
                    suppressWarnings() |>
                    suppressMessages()

                ests[, y, r] = exp(predict(m, type="link"))
            }

            # normalize
            ests[, ref_grp, r] = 1
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
        ests_list = lapply(seq_len(n_grp), function(j) ests[j, , ])
        calc_bayes(Y, idx_grp, ests_list, p_rxs, n_grp)

        if (break_next) break
        if (all(!updating)) { # converged -- do one more full update
            updating = matrix(TRUE, n_y, n_r)
            break_next = TRUE
        }
    }
    cli::cli_progress_done()


    # final global mean
    est = dirichlet_map(Y, ones, p_ryxs, rep_len(1, n_y), 1)[[1]]

    list(map = est,
         ests = ests,
         max_iters = i >= iter,
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

