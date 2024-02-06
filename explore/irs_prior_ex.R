library(birdie)
# `morg_int` is the outcome variable
# `r_probs` contains the BISG estiamtes
# `itm_data` is the data frame

# empirical bayes prior
# 0.1 is the effective number of data points per geo unit (PUMA/ZCTA) in the prior:
# make it bigger for the prior to have a stronger regularizing effect
prior_eb = list(
    alpha = 1 + 0.1 * coef(est_weighted(r_probs, morg_int ~ 1, data=itm_data))
)

# then use e.g. as follows:
fit1 = birdie(r_prob, morg_int ~ puma, data=itm_data, prior=prior_eb)

# passing boot = 10 to the above will do the bootstrapping we talked about as well
