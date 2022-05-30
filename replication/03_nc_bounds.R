p_xr_party_true = prop.table(table(d$party, d$race))
p_xr_turnout_true = prop.table(table(d$n_voted, d$race))
p_xr_party_est = calc_joint_model(fits_party$county, q=0.5, p_r=p_r)
p_xr_turnout_est = calc_joint_model(fits_turnout$county, q=0.5, p_r=p_r)

p_xr_party_low = calc_joint_model(fits_party$block, q=0.05, p_r=p_r)
p_xr_party_high = calc_joint_model(fits_party$block, q=0.95, p_r=p_r)
idx_max = which.max(p_xr_party_high - p_xr_party_low)

list(
    ci_min = p_xr_party_low[idx_max],
    ci_max = p_xr_party_high[idx_max],
    tot_err_party = max(abs(p_xr_party_true - p_xr_party_est) /
                            fits_party$county$raw_error$global),
    tot_err_turnout = max(abs(p_xr_turnout_true - p_xr_turnout_est) /
                              fits_turnout$county$raw_error$global),

    avg_err_party = tot_err_party / sqrt(nrow(d)),
    avg_err_turnout = tot_err_turnout / sqrt(nrow(d))
) |>
    write_rds(here("paper/data/nc_bounds.rds"))

