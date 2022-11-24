suppressMessages(devtools::load_all(here::here(".")))

data(pseudo_vf)
pseudo_vf = do.call(rbind, lapply(1:100, \(i) pseudo_vf))

r_probs = bisg(~ nm(last_name), pseudo_vf) |>
    as.matrix()

Y = as.integer(pseudo_vf$turnout)

# ones = rep_along(Y, 1L)
X = to_unique_ids(coalesce(pseudo_vf$zip, "none"))
n_x = max(X)

est0 = dirichlet_map(Y, X, r_probs, c(1, 1), n_x)

step2 = function() {
    pr = calc_bayes(Y, X, est0, r_probs, n_x, 2)
    dirichlet_map(Y, X, pr, c(1, 1), n_x)
}


all.equal(step2(), em_dirichlet(est0, Y, X, r_probs, c(1, 1), n_x))

microbenchmark::microbenchmark(
    step2(),
    em_dirichlet(est0, Y, X, r_probs, c(1, 1), n_x)
)


# x = birdie(r_probs, turnout ~ 1, pseudo_vf)
# print(x)

# x = birdie(r_probs, turnout ~ (1 | zip), pseudo_vf, max_iter=50)
# print(x)
