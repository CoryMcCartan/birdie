

p_xr = prop.table(table(pseudo_vf$turnout, pseudo_vf$race), 2)
p_r1 = prop.table(table(pseudo_vf$race))
p_r2 = p_r1 * c(1.15, 0.95, 0.9, 0.9, 3.0, 0.8)
p_r2 = p_r2 / sum(p_r2)
p_r2-p_r1

r_pr1 = bisg(~ nm(last_name) + zip(zip), data=pseudo_vf, p_r=p_r1) |> # true
    as.matrix()
r_pr2 = bisg(~ nm(last_name) + zip(zip), data=pseudo_vf, p_r=p_r2) |> # wrong
    as.matrix()

t(r_pr2) %*% (r_pr2 - r_pr1) %*% p_xr[2, ]
solve(t(r_pr2) %*% r_pr2) %*% t(r_pr2) %*% (r_pr1 - r_pr2) %*% p_xr[2, ]

coef(lm(pseudo_vf$turnout == "yes" ~ r_pr2 + 0)) - coef(lm(pseudo_vf$turnout == "yes" ~ r_pr1 + 0))
coef(lm(pseudo_vf$turnout == "yes" ~ r_pr2 + 0)) - p_xr[2,]


r_probs = bisg(~ nm(last_name) + zip(zip), pseudo_vf)
p_r = prop.table(table(pseudo_vf$race))
p_r_est = colMeans(r_probs)
X = as.matrix(r_probs)
mu = lm.fit(X, pseudo_vf$turnout == "yes")$coefficients

invM = solve(t(X) %*% X)
sens_r = map_dbl(1:6, function(i) {
    delta = matrix(0, nrow=nrow(X), ncol=ncol(X))
    delta[, i] = 0.5 * X[, i]
    delta_wt = 0.9*X[, -i] + 0.1 * 1/5
    delta[, -i] = -delta[, i] * delta_wt / rowSums(delta_wt)
    (invM %*% t(X) %*% delta %*% mu)[i, 1]
})
abs(sens_r / sens_r[1])
