suppressMessages(devtools::load_all(here::here(".")))

data(pseudo_vf)

r_probs = bisg(~ nm(last_name), pseudo_vf)

Y = as.integer(pseudo_vf$turnout)

x = birdie(r_probs, turnout ~ 1, pseudo_vf)
print(x)

# x = birdie(r_probs, turnout ~ (1 | zip), pseudo_vf, max_iter=50)
# print(x)
