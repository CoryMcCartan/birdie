p_xr_party_low = calc_joint_model(fits_party$zip, q=0.05, p_r=p_r)
p_xr_party_high = calc_joint_model(fits_party$zip, q=0.95, p_r=p_r)
idx_max = which.max(p_xr_party_high - p_xr_party_low)

err_wb = sqrt(6 * nrow(d) * rowMeans((fits_party$county$raw_cov$global[1, , ]
                                      - fits_party$county$raw_cov$global[2, , ])^2))
err_wh = sqrt(6 * nrow(d) * rowMeans((fits_party$county$raw_cov$global[1, , ]
                                      - fits_party$county$raw_cov$global[3, , ])^2))

disp_wb = xr_party$model_zip %*% diag(1/p_r) %*% c(1, -1, 0, 0, 0, 0)
disp_wh = xr_party$model_zip %*% diag(1/p_r) %*% c(1, 0, -1, 0, 0, 0)

Delta = seq(0, 12, by=0.5)
total_to_avg = 1 / sqrt(nrow(d))

d_sens_wb = map_dfr(Delta, function(x) {
    cbind(disp_wb, err_wb) |>
        as.data.frame() |>
        rownames_to_column("party") |>
        rename(est=V1) |>
        mutate(delta_tot = x,
               delta_avg = x * total_to_avg,
               party = str_to_upper(party),
               low = est - x * err_wb,
               high = est + x * err_wb)
})
d_sens_wh = map_dfr(Delta, function(x) {
    cbind(disp_wh, err_wh) |>
        as.data.frame() |>
        rownames_to_column("party") |>
        rename(est=V1) |>
        mutate(delta_tot = x,
               delta_avg = x * total_to_avg,
               party = str_to_upper(party),
               low = est - x * err_wb,
               high = est + x * err_wb)
})

p1 = ggplot(d_sens_wb, aes(delta_avg, est, color=party, fill=party)) +
    geom_line(linetype="dashed") +
    geom_ribbon(aes(ymin=low, ymax=high), alpha=0.5, color=NA) +
    geom_hline(yintercept=0) +
    scale_x_continuous("Average individual BISG error, ||δ||",
                       labels=label_percent(0.1, suffix="pp")) +
    scale_y_continuous("Estimated disparity",
                       labels=label_percent(1.0, suffix="pp")) +
    scale_fill_manual(values=c(LIB="#E4C22B", REP="#B25D4C", IND="#aaccaf", DEM="#3D77BB")) +
    scale_color_manual(values=c(LIB="#E4C22B", REP="#B25D4C", IND="#aaccaf", DEM="#3D77BB")) +
    coord_cartesian(ylim=c(-1, 1), expand=FALSE) +
    labs(title="White-Black disparity", color="Party", fill="Party") +
    theme_paper()
p2 = ggplot(d_sens_wh, aes(delta_avg, est, color=party, fill=party)) +
    geom_line(linetype="dashed") +
    geom_ribbon(aes(ymin=low, ymax=high), alpha=0.5, color=NA) +
    geom_hline(yintercept=0) +
    scale_x_continuous("Average individual BISG error, ||δ||",
                       labels=label_percent(0.1, suffix="pp")) +
    scale_y_continuous("Estimated disparity",
                       labels=label_percent(1.0, suffix="pp")) +
    scale_fill_manual(values=c(LIB="#E4C22B", REP="#B25D4C", IND="#aaccaf", DEM="#3D77BB")) +
    scale_color_manual(values=c(LIB="#E4C22B", REP="#B25D4C", IND="#aaccaf", DEM="#3D77BB")) +
    coord_cartesian(ylim=c(-1, 1), expand=FALSE) +
    labs(title="White-Hispanic disparity", color="Party", fill="Party") +
    theme_paper()

p = p1 + p2 + plot_layout(guides="collect")
ggsave(here("paper/figures/nc_bounds.pdf"), plot=p, width=7.5, height=4, device=cairo_pdf)


list(
    ci_min = p_xr_party_low[idx_max],
    ci_max = p_xr_party_high[idx_max],
    max_err_party_wb = abs(disp_wb/err_wb)[,1] * total_to_avg,
    max_err_party_wh = abs(disp_wh/err_wh)[,1] * total_to_avg,
    total_to_avg = total_to_avg
) |>
    write_rds(here("paper/data/nc_bounds.rds"))

