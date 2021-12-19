library(raceproxy)

# setup ------
# a data frame containing address info and the outcome variable
tax_data = ...

# should be a data frame with a row for each individual
# and columns pr_white pr_black pr_hisp pr_asian pr_aian pr_other
r_probs = ...

# these are our best guesses for the total racial distribution in the sample
# we can pull this from the CPS estimate of the taxpaying population,
# if the IRS doesn't have its own estimates
p_r = c(white=0.630, black=0.121, hisp=0.173,
        asian=0.0478, aian=0.0072, other=0.0210)

# run the procedure -----
# start with this. If the `rhat` value in the output is getting stuck & not
# going lower, may try increasing `lr` to e.g. 0.3 or 0.4.
config = list(it_avg=400, tol_rhat=1.1, lr=0.25)
# replace `X` with the outcome variable
# replace `ADDRESS` with the address variable (zip, Census block, etc.)
fit = model_race(r_probs, X, ADDRESS, data=tax_data, config=config)

# check some diagnostics ----
# these should look like they converge to some asymptote
plot(fit$loss, type='l')
plot(fit$loss[-1:-20], type='l')

# model output -------
# will be number of approximate posterior draws BY number of outcome levels
# BY number of racial categories
dim(fit$p_xr)

# the point estimate for the outcome-by-race joint table
est_xr = calc_joint_model(fit$p_xr, 0.5, p_r)

# credible intervals. These will be very anticonservative in practice because of
# data biases.
est_xr_low = calc_joint_model(fit$p_xr, 0.05, p_r)
est_xr_high = calc_joint_model(fit$p_xr, 0.95, p_r)

# convert from joint table to make each column (race group) sum to 1
est_xr %*% diag(1/p_r)
