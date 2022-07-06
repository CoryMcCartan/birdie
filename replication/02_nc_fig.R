eval_fit = function(fits, X) {
    X = eval_tidy(enquo(X), d)

    xr = list(true = prop.table(table(X, d$race)))
    xr = c(xr, flatten(imap(fits, function(f, level) {
        out = list()
        out[[str_c("model_", level)]] = calc_joint_model(f, "global", 0.5, p_r)
        out
    })))
    xr = c(xr, flatten(imap(r_probs, function(d_pr, level) {
        out = list()
        out[[str_c("weight_", level)]] = calc_joint_bisgz(d_pr, X, method="weight")
        out[[str_c("thresh_", level)]] = calc_joint_bisgz(d_pr, X, method="thresh")
        out[[str_c("ols_", level)]] = calc_joint_bisgz(d_pr, X, method="ols")
        out
    })))

    tv_overall = do.call(eval_joints, c(list(xr$true, "tv"), xr)) %>%
        filter(method != "true") %>%
        mutate(race = "overall") %>%
        rename(tv=TV)

    tv_race = do.call(eval_joints, c(list(xr$true, "tv_col"), xr)) %>%
        filter(method != "true") %>%
        unnest_longer(TV_COL) %>%
        rename(tv=TV_COL, race=TV_COL_id)

    bind_rows(tv_overall, tv_race) %>%
        separate(method, c("method", "level"), sep="_")
}

tv_party = eval_fit(fits_party, party)
tv_turnout = eval_fit(fits_turnout, n_voted)

# Overview plots -----

p1 = ggplot(d, aes(y=factor(races[race], levels=rev(races)),
              fill=factor(str_to_upper(party), levels=c("DEM", "IND", "REP", "LIB")))) +
    geom_bar(position="fill") +
    scale_fill_manual(values=c(LIB="#E4C22B", REP="#B25D4C", IND="#aaccaf", DEM="#3D77BB")) +
    scale_x_continuous("Proportion", labels=percent, expand=c(0, 0)) +
    labs(y=NULL, fill="Party") +
    theme_paper()

p2 = ggplot(d, aes(y=factor(races[race], levels=rev(races)),
                   fill=as.numeric(as.character(n_voted)),
                   group=fct_rev(n_voted))) +
    geom_bar(position="fill") +
    scale_fill_wa_c("forest_fire", breaks=c(1, 3, 5, 7, 9)) +
    scale_x_continuous("Proportion", labels=percent, expand=c(0, 0)) +
    labs(y=NULL, fill="Elections voted,\n2012-2021") +
    theme_paper()

p = p1 + p2 & theme(legend.position="bottom",
                    legend.margin=margin())
ggsave(here("paper/figures/nc_overview.pdf"), plot=p, width=7.5, height=3.5)

# Fit quality plots -----

lbl_race = function(race) {
    str_c(races[race], " (", percent(as.numeric(p_r[race]), 0.1), ")")
}

geos = c(county="County", zip="ZIP code", tract="Census tract", block="Census block")
geos_short = c(county="County", zip="ZIP", tract="Tract", block="Block")
methods = c(model="Model", ols="OLS", thresh="Threshold", weight="Weighted")

## Party fit quality -----
p1 = filter(tv_party, race=="overall") %>%
ggplot(aes(x=factor(geos[level], levels=geos), y=tv,
           color=methods[method], shape=methods[method], group=methods[method])) +
    geom_textline(aes(label=methods[method]),
                  position=position_dodge(width=0.25),
                  linewidth=0.7, size=3.5, family="Times", hjust=0.08) +
    geom_point(size=2.0, position=position_dodge(width=0.25)) +
    scale_color_wa_d() +
    scale_x_discrete(expand=c(0.07, 0, 0.06, 0)) +
    scale_y_log10("Overall total variation distance") +
    labs(x="BISG geographic precision", title="Party identification") +
    guides(color="none", shape="none") +
    theme_paper() +
    theme(plot.margin=unit(c(0, 0.1, 0, 0), "cm"),
          plot.title=element_text(margin=margin(0, 0, -12, 0)))

p2 = filter(tv_party, race!="overall") %>%
ggplot(aes(x=factor(geos_short[level], levels=geos_short), y=tv,
           color=methods[method], shape=methods[method], group=methods[method])) +
    facet_wrap(~ factor(lbl_race(race), levels=lbl_race(names(races)))) +
    geom_line(position=position_dodge(width=0.25), size=0.7) +
    geom_point(size=2.0, position=position_dodge(width=0.25)) +
    scale_color_wa_d() +
    scale_y_log10("Total variation distance") +
    labs(x="BISG geographic precision") +
    guides(color="none", shape="none") +
    theme_paper() +
    theme(plot.margin=unit(c(0, 0, 0, 0), "cm"))

p = p1 + p2 + plot_layout(widths=c(0.4, 0.6))
ggsave(here("paper/figures/nc_party_fit.pdf"), plot=p, width=8, height=3.75)

## Turnout fit quality -----
p1 = filter(tv_turnout, race=="overall") %>%
ggplot(aes(x=factor(geos[level], levels=geos), y=tv,
           color=methods[method], shape=methods[method], group=methods[method])) +
    geom_textline(aes(label=methods[method],
                      hjust=c(model=0.08, weight=0.08, thresh=0.52, ols=0.1)[method]),
                  position=position_dodge(width=0.25),
                  linewidth=0.7, size=3.5, family="Times") +
    geom_point(size=2.0, position=position_dodge(width=0.25)) +
    scale_color_wa_d() +
    scale_x_discrete(expand=c(0.07, 0, 0.06, 0)) +
    scale_y_log10("Overall total variation distance") +
    labs(x="BISG geographic precision", title="Turnout") +
    guides(color="none", shape="none") +
    theme_paper() +
    theme(plot.margin=unit(c(0, 0.1, 0, 0), "cm"),
          plot.title=element_text(margin=margin(0, 0, -12, 0)))

p2 = filter(tv_turnout, race!="overall") %>%
ggplot(aes(x=factor(geos_short[level], levels=geos_short), y=tv,
           color=methods[method], shape=methods[method], group=methods[method])) +
    facet_wrap(~ factor(lbl_race(race), levels=lbl_race(names(races)))) +
    geom_line(position=position_dodge(width=0.25), size=0.7) +
    geom_point(size=2.0, position=position_dodge(width=0.25)) +
    scale_color_wa_d() +
    scale_y_log10("Total variation distance") +
    labs(x="BISG geographic precision") +
    guides(color="none", shape="none") +
    theme_paper() +
    theme(plot.margin=unit(c(0, 0, 0, 0), "cm"))

p = p1 + p2 + plot_layout(widths=c(0.4, 0.6))
ggsave(here("paper/figures/nc_turnout_fit.pdf"), plot=p, width=8, height=3.75)


# Poster plots -----
# bind_rows(`Party ID`=tv_party, `Turnout`=tv_turnout, .id="outcome") |>
tv_party |>
    filter(method %in% c("model", "weight"), level == "zip") |>
    mutate(race = fct_inorder(c(overall="Overall", races)[race])) |>
ggplot(aes(race, tv, fill=c(model="New model", weight="Weighting")[method])) +
    geom_hline(yintercept=0) +
    geom_col(position=position_dodge()) +
    labs(title="Party registration", x=NULL,
         y="TV distance\n(lower is better)", fill="Method") +
    scale_fill_wa_d("palouse", which=c("snake", "hills")) +
    coord_cartesian(expand=FALSE) +
    theme_minimal(base_family="IBM Plex Sans Medium", base_size=24) +
    theme(plot.margin=margin(0, 0, 0, 0),
          legend.position=c(0.8, 0.55),
          panel.grid.major.x=element_blank())
ggsave("~/Desktop/chart1.png", width=12, height=4.5, dpi=250)

tv_party |>
    filter(method %in% c("model", "weight"), level == "zip") |>
    mutate(race = fct_inorder(c(overall="Overall", races)[race])) |>
ggplot(aes(race, tv, fill=c(model="New model", weight="Weighting")[method])) +
    geom_hline(yintercept=0) +
    geom_col(position=position_dodge()) +
    labs(title="Turnout", x=NULL, y=NULL, fill="Method") +
    scale_fill_wa_d("palouse", which=c("snake", "hills"), guide="none") +
    coord_cartesian(expand=FALSE) +
    theme_minimal(base_family="IBM Plex Sans Medium", base_size=24) +
    theme(plot.margin=margin(0, 0, 0, 0),
          panel.grid.major.x=element_blank())
ggsave("~/Desktop/chart2.png", width=11, height=4.5, dpi=250)
