make_est_d = function(fits, X, d) {
    X_name = deparse(substitute(X))
    X = eval_tidy(enquo(X), d)

    d_true = tidy_true(d$race, X, X_name)

    d_birdie = map(fits, function(ll) {
        bind_rows(map(ll, tidy), .id = "method")
    }) |>
        bind_rows(.id = "level") |>
        mutate(method = str_c("birdie_", method))

    form = as.formula(paste(X_name, "~ 1"))
    d_wtd = map(r_probs, function(d_pr) {
        tidy(est_weighted(d_pr, form, data=d))
    }) |>
        bind_rows(.id = "level") |>
        mutate(method = "weight")

    # d_ols = map(r_probs, function(d_pr) {
    #     tidy_ols(d_pr, X, X_name)
    # }) |>
    #     bind_rows(.id = "level") |>
    #     mutate(method = "ols")

    d_thresh = map(r_probs, function(d_pr) {
        tidy_thresh(d_pr, X, X_name)
    }) |>
        bind_rows(.id = "level") |>
        mutate(method = "thresh")

    bind_rows(d_birdie, d_wtd,  d_thresh) |>
        rename(est = estimate) |>
        left_join(d_true, by=c(X_name, "race")) |>
        rename(est_true = estimate)
}

calc_disp <- function(ests) {
    ests |>
        filter(race %in% c("white", "black", "hisp")) |>
        pivot_wider(names_from=race, values_from=c(est, est_true)) |>
        mutate(disp_wb = est_white - est_black,
               disp_wb_true = est_true_white - est_true_black,
               disp_wh = est_white - est_hisp,
               disp_wh_true = est_true_white - est_true_hisp) |>
        select(-starts_with("est_"))
}

eval_fit_tv = function(ests, p_r) {
    tv_overall = ests |>
        group_by(level, method) |>
        summarize(tv = sum(abs(est - est_true)*p_r[race])/2) |>
        mutate(race = "overall")

    tv_race = ests |>
        group_by(level, method, race) |>
        summarize(tv = sum(abs(est - est_true))/2)

    ungroup(bind_rows(tv_overall, tv_race))
}


if (!file.exists(path <- here("data-out/nc_ests_party.rda"))) {
    ests_party = make_est_d(fits_party, party, d)
    save(ests_party, p_r, file=path, compress="gzip")
    write_rds(ests_party, path, compress="gz")
} else {
    load(path)
    ests_party = read_rds(path)
}

disp_party = calc_disp(ests_party)

tv_party = eval_fit_tv(ests_party, p_r)

filter(disp_party, level=="zip", party=="dem") |>
    select(-level, -party) |>
    split(~ method) |>
    write_rds(here("paper/data/nc_disp_ex.rds"), compress="xz")

# Overview plots -----

p1 = ggplot(d, aes(y=factor(races[race], levels=rev(races)),
              fill=factor(str_to_upper(party), levels=c("DEM", "IND", "REP", "LIB")))) +
    geom_bar(position="fill") +
    scale_fill_manual(values=c(LIB="#E4C22B", REP="#B25D4C", IND="#aaccaf", DEM="#3D77BB")) +
    scale_x_continuous("Proportion", labels=percent, expand=c(0, 0)) +
    labs(y=NULL, fill="Party") +
    theme_paper()


ggsave(here("paper/figures/nc_overview.pdf"), plot=p1, width=6.5, height=2.5)

# Disparity plots -----

lbl_race = function(race) {
    str_c(races[race], " (", percent(as.numeric(p_r[race]), 0.1), ")")
}

geos = c(county="County", zip="ZIP code", tract="Census tract", block="Census block")
geos_short = c(county="County", zip="ZIP", tract="Tract", block="Block")
methods = c(birdie_pool="BIRDiE (pooling)", birdie_sat="BIRDiE (saturated)",
            birdie_mmm="BIRDiE (mixed)",
            ols="OLS", thresh="Threshold", weight="Weighting")
methods_short = c(birdie_pool="BIRDiE (pool.)", birdie_sat="BIRDiE (sat.)",
                  birdie_mmm="BIRDiE (mixed)",
                  ols="OLS", thresh="Threshold", weight="Weighting")
methods_col = c(birdie_pool=PAL_R[3], birdie_sat=PAL_R[3], birdie_mmm=PAL_R[3],
                ols=PAL_R[2], thresh=PAL_R[4], weight=PAL_R[1])
methods_shp = c(birdie_pool=16, birdie_sat=15, birdie_mmm=16,
                ols=3, thresh=4, weight=1)


d_ann <- tibble(
    x = 1.4,
    y = with(disp_party, disp_wb_true[party == "dem"][1]),
    label = "True disparity",
    groups = "wb"
)
filter(disp_party, level == "block") |>
    rename(disp_wb_est=disp_wb, disp_wh_est=disp_wh) |>
    pivot_longer(starts_with("disp"), names_pattern="disp_(w.)_(.+)", names_to=c("groups", "version")) |>
    pivot_wider(names_from=version) |>
ggplot(aes(str_to_upper(party), est, color=method, shape=method, group=method)) +
    facet_wrap(~ groups, labeller = \(...) list(
        groups=c("White-Black disparity", "White-Hispanic disparity")
    )) +
    geom_hline(yintercept=0.0, color="#00000077") +
    geom_blank() +
    geom_segment(aes(x=as.integer(as.factor(party))-0.35,
                     xend=as.integer(as.factor(party))+0.35,
                     y=true, yend=true),
                 lty="11", col="#444444", linewidth=0.6) +
    geom_point(size=2.4, position=position_dodge(width=0.7)) +
    geom_text(aes(x=x, y=y, label=label), data=d_ann, inherit.aes=FALSE,
              hjust=0, size=2.5, family="Times") +
    scale_y_continuous("Disparity estimate",
                       labels=label_number(1, scale=100, suffix="pp")) +
    scale_color_manual(values=methods_col, labels=methods) +
    scale_shape_manual(values=methods_shp, labels=methods) +
    labs(x="Party", color="Method", shape="Method") +
    theme_paper()

ggsave(here("paper/figures/nc_disp.pdf"), width=8, height=3.5) # old h 5



# TV plots -----

## Party fit quality -----
p1 = filter(tv_party, race=="overall") %>%
ggplot(aes(x=factor(geos[level], levels=geos), y=tv,
           color=method, shape=method, group=method)) +
    geom_textline(aes(label=methods_short[method],
                      hjust=c(birdie_sat=0.08, birdie_pool=0.08, birdie_mmm=0.50,
                              weight=0.08, thresh=0.08)[method]),
                  position=position_dodge(width=0.7),
                  linewidth=0.7, size=3.0, family="Times") +
    geom_point(size=2.4, position=position_dodge(width=0.7)) +
    scale_x_discrete(expand=c(0.07, 0, 0.06, 0)) +
    scale_y_continuous("Overall total variation distance",
                       limits=c(0.0, 0.13), expand=expansion(c(0, 0.05))) +
    # labs(x="BISG geographic precision", title="Party Identification") +
    labs(x="BISG geographic precision") +
    # scale_color_manual(values=methods_col, labels=methods, guide="none") +
    # scale_shape_manual(values=methods_shp, labels=methods, guide="none") +
    scale_color_manual(values=methods_col, labels=methods) +
    scale_shape_manual(values=methods_shp, labels=methods) +
    labs(color=NULL, shape=NULL) +
    guides(color=guide_legend(override.aes=aes(label = ""))) +
    theme_paper() +
    theme(plot.margin=unit(c(0, 0.1, 0, 0), "cm"),
          legend.text=element_text(size=8.0),
          legend.key.height=unit(0.4, "cm"),
          legend.background=element_blank(),
          legend.position=c(0.7, 0.37))

p1b = filter(tv_party, race!="overall") %>%
ggplot(aes(x=factor(geos_short[level], levels=geos_short), y=tv,
           color=method, shape=method, group=method)) +
    facet_wrap(~ factor(lbl_race(race), levels=lbl_race(names(races)))) +
    geom_line(position=position_dodge(width=0.7), linewidth=0.7) +
    geom_point(size=2.0, position=position_dodge(width=0.7)) +
    scale_y_continuous("Total variation distance", limits=c(0, 0.31),
                       expand=expansion(c(0, 0.03))) +
    labs(x="BISG geographic precision", #title="Party",
         color="Method", shape="Method") +
    scale_color_manual(values=methods_col, labels=methods) +
    scale_shape_manual(values=methods_shp, labels=methods) +
    guides(shape="none", color="none") +
    theme_paper() +
    theme(plot.margin=unit(c(0, 0, 0, 0), "cm"))


p = p1 + p1b
ggsave(here("paper/figures/nc_tv.pdf"), plot=p, width=8, height=4) # old w 8

# p = p1 + p1b
# ggsave(here("paper/figures/nc_tv_detailed.pdf"), plot=p, width=8, height=5) # old h 9


# Small-area estimates -----------

d_small_true = map_dfr(geo_levels[1:3], function(l) {
    d |>
        rename(GEOID=str_c("GEOID_", l)) |>
        count(party, race, GEOID) |>
        group_by(race, GEOID) |>
        mutate(est_true = n / sum(n)) |>
        rename(pop = n) |>
        ungroup() |>
        mutate(level = l, .before=everything())
})

d_small = map_dfr(geo_levels[1:3], function(l) {
    d_lv = rename(d, GEOID=str_c("GEOID_", l))
    idxs = split(seq_len(nrow(d)), d_lv$GEOID)
    bind_rows(
        birdie_mmm = tidy(fits_party[[l]]$mmm, subgroup=TRUE),
        birdie_sat = tidy(fits_party[[l]]$sat, subgroup=TRUE),
        weight = tidy(
            est_weighted(r_probs[[l]], party ~ GEOID, data=d_lv),
            subgroup=TRUE),
        thresh = map(idxs, ~ tidy_thresh(r_probs[[l]][., ], d$party[.], "party")) |>
            bind_rows(.id="GEOID"),
        .id="method"
    ) |>
        mutate(level = l, .before=everything())
}) |>
    select(-any_of(c("white", "black"))) |>
    rename(est=estimate) |>
    left_join(d_small_true, by=c("party", "race", "level", "GEOID")) |>
    filter(pop >= 5) |>
    drop_na(est) |>
    mutate(race = fct_inorder(race),
           level = fct_inorder(geos_short[level]),
           est_true = coalesce(est_true, 0)) |>
    group_by(method, level, race, GEOID) |>
    mutate(tv = sum(abs(est - est_true)) / 2,
           rmse = sqrt(mean((est - est_true)^2)),
           pop = sum(pop)) |>
    suppressWarnings() |>
    group_by(method, level, race) |>
    summarize(tv_wt = weighted.mean(tv, pop),
              tv = mean(tv),
              rmse = mean(rmse),
              cor_all = cor(est, est_true, method="spearman"),
              .groups="drop")

p1 = d_small |>
    filter(race %in% c("white", "black")) |>
ggplot(aes(level, tv, color=method, shape=method, group=method)) +
    facet_grid(~ fct_inorder(lbl_race(race))) +
    geom_line(linewidth=0.7, position=position_dodge(0.7)) +
    geom_point(size=2.0, position=position_dodge(0.7)) +
    scale_color_manual(values=methods_col, labels=methods_short) +
    scale_shape_manual(values=methods_shp, labels=methods_short) +
    scale_y_continuous(limits=c(0, NA), expand=expansion(mult=c(0, 0.05))) +
    labs(x="BISG geographic precision", y="Mean TV distance across areas", color="Method", shape="Method") +
    theme_paper() +
    theme(legend.margin=margin(),
          plot.margin=margin())

p2 = d_small |>
    filter(race %in% c("white", "black")) |>
ggplot(aes(level, rmse, color=method, shape=method, group=method)) +
    facet_grid(~ fct_inorder(lbl_race(race))) +
    geom_line(linewidth=0.7, position=position_dodge(0.7)) +
    geom_point(size=2.0, position=position_dodge(0.7)) +
    scale_color_manual(values=methods_col, labels=methods_short) +
    scale_shape_manual(values=methods_shp, labels=methods_short) +
    scale_y_continuous(limits=c(0, NA), expand=expansion(mult=c(0, 0.05))) +
    labs(x="BISG geographic precision", y="Mean RMSE across areas", color="Method", shape="Method") +
    theme_paper() +
    theme(legend.margin=margin(),
          plot.margin=margin())
p3 = d_small |>
    filter(race %in% c("white", "black")) |>
ggplot(aes(level, cor_all, color=method, shape=method, group=method)) +
    facet_grid(~ fct_inorder(lbl_race(race))) +
    geom_line(linewidth=0.7, position=position_dodge(0.7)) +
    geom_point(size=2.0, position=position_dodge(0.7)) +
    scale_color_manual(values=methods_col, labels=methods_short) +
    scale_shape_manual(values=methods_shp, labels=methods_short) +
    labs(x="BISG geographic precision", y="Correlation with truth, all estimates",
         color="Method", shape="Method") +
    theme_paper() +
    theme(legend.margin=margin(),
          plot.margin=margin(l=12))

ggsave(here("paper/figures/nc_smallarea.pdf"), plot=p1, width=6.5, height=2.75)

p = p2 + p3 + plot_layout(nrow=1, guides="collect")
ggsave(here("paper/figures/nc_smallarea_app.pdf"), plot=p, width=8, height=3)
