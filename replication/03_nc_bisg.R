# Updated BISG quality ----------
log_scores_party = map_dbl(fits_party$block, function(x) {
    x = as.matrix(fitted(x))
    pr_act = x[cbind(1:nrow(x), as.integer(d$race))]
    pr_act[pr_act == 0] = 1e-6
    mean(log(pr_act))
})
acc_thresh_party = map_dbl(map(fits_party$block, fitted),
                           ~ mean(max.col(.) == as.integer(d$race)))

list(score = log_scores_party,
     acc = acc_thresh_party) |>
    write_rds(here("paper/data/nc_bisg_post_party.rds"))

d_roc = map_dfr(geo_levels, function(l) {
    cat("Calculating AUROC at the", l, "level\n")
    bind_rows(
        bisg = calc_roc(r_probs[[l]]),
        sat = calc_roc(fitted(fits_party[[l]]$sat)),
        mmm = calc_roc(fitted(fits_party[[l]]$mmm)),
        .id="method"
    ) |>
        mutate(level = l, .before=everything())
})

methods = c(sat="BIRDiE (sat.)", mmm="BIRDiE (mixed)", bisg="BISG")
methods_col = c(sat=PAL_R[3], mmm=PAL_R[3], bisg=PAL_R[2])
methods_shp = c(sat=15, mmm=16, bisg=17)

d_roc |>
    pivot_longer(white:other, names_to="race", values_to="roc") |>
ggplot(aes(fct_inorder(geos_short[level]), roc,
           color=method, shape=method, group=method)) +
    facet_grid(~ fct_inorder(lbl_race(race))) +
    geom_line(linewidth=0.7, position=position_dodge(0.7)) +
    geom_point(size=2.0, position=position_dodge(0.7)) +
    scale_color_manual(values=methods_col, labels=methods) +
    scale_shape_manual(values=methods_shp, labels=methods) +
    scale_y_continuous("Area under the ROC curve") +
    labs(x="Geographic precision", color="Method", shape="Method") +
    theme_paper() +
    theme(axis.text.x=element_text(size=7),
          plot.margin=margin(),
          legend.margin=margin())

ggsave(here("paper/figures/nc_roc.pdf"), width=8.0, height=2.25)

