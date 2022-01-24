est_bisg = function(X, S, G, p_sr, p_gzr, p_r, geo=TRUE) {
    if (!geo) p_gzr = p_gzr*0 + 1

    m_bisg = matrix(nrow=length(X), ncol=length(p_r))
    for (i in seq_along(p_r)) {
        m_bisg[, i] = p_sr[S, i] * p_gzr[G, i] * p_r[i]
    }
    m_bisg = m_bisg / rowSums(m_bisg)
    colnames(m_bisg) = names(p_r)
    m_bisg
}

est_leastsq = function(X, S, p_sr, p_r) {
    rlang::check_installed(c("MASS", "glmnet"), "generalized inverse")
    p_rs = p_sr %*% diag(p_r)
    p_rs_inv = t(MASS::ginv(p_rs))
    p_xs = prop.table(table(X, S))

    p_x = prop.table(table(X))
    P_xr_low = outer(p_x, p_r, \(x, y) pmax(0, x+y-1))
    P_xr_high = outer(p_x, p_r, pmin)

    tidy_pred = function(est) {
        est[est > P_xr_high] = P_xr_high[est > P_xr_high]
        est[est < P_xr_low] = P_xr_low[est < P_xr_low]
        rownames(est) = levels(X)
        colnames(est) = names(p_r)
        est
    }

    P_xr_lsq = tidy_pred(p_xs %*% p_rs_inv)
    P_xr_nnls = tidy_pred(rbind(
        as.numeric(glmnet::glmnet(p_rs, t(p_xs[1, ]), lambda=0, lower.limits=0, intercept=FALSE)$beta),
        as.numeric(glmnet::glmnet(p_rs, t(p_xs[2, ]), lambda=0, lower.limits=0, intercept=FALSE)$beta),
        as.numeric(glmnet::glmnet(p_rs, t(p_xs[3, ]), lambda=0, lower.limits=0, intercept=FALSE)$beta),
        as.numeric(glmnet::glmnet(p_rs, t(p_xs[4, ]), lambda=0, lower.limits=0, intercept=FALSE)$beta)
    ) %*% diag(p_r))

    list(lsq = P_xr_lsq,
         nnls = P_xr_nnls)
}
