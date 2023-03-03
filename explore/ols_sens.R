

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
