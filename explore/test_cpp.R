suppressMessages(devtools::load_all(here::here(".")))

data(pseudo_vf)
# pseudo_vf = do.call(rbind, lapply(1:100, \(i) pseudo_vf))

r_probs = bisg(~ nm(last_name), pseudo_vf) |>
    as.matrix()

Y = as.integer(pseudo_vf$turnout)

ones = rep_along(Y, 1L)
est0 = dirichlet_map(Y, ones, r_probs, c(1, 1), 1)

str(est0)

# x = birdie(r_probs, turnout ~ 1, pseudo_vf)
# print(x)

# x = birdie(r_probs, turnout ~ (1 | zip), pseudo_vf, max_iter=50)
# print(x)
