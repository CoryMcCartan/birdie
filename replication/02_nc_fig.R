eval_fit_disp = function(fits, X) {
    X_name = deparse(substitute(X))
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

    m_to_cond = diag(1/p_r)
    rownames(m_to_cond) = names(p_r)
    colnames(m_to_cond) = names(p_r)

    x_disp = imap_dfr(xr, function(x, nm) {
        x = x %*% m_to_cond
        tibble(alg = nm,
               {{ X_name }} := rownames(x),
               disp_wb = x[, "white"] - x[, "black"],
               disp_wh = x[, "white"] - x[, "hisp"])
    })
    x_disp_true = filter(x_disp, alg == "true") |>
        select(-alg)
    x_disp |>
        filter(alg != "true") |>
        separate(alg, c("method", "level"), sep="_") |>
        left_join(x_disp_true, by=X_name, suffix=c("", "_true"))
}

eval_fit_tv = function(fits, X) {
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

disp_party = eval_fit_disp(fits_party, party)
disp_turnout = eval_fit_disp(fits_turnout, n_voted)

tv_party = eval_fit_tv(fits_party, party)
tv_turnout = eval_fit_tv(fits_turnout, n_voted)

bind_rows(party=tv_party, turnout=tv_turnout, .id="outcome") |>
    filter(race == "overall", method %in% c("model", "weight")) |>
    pivot_wider(names_from=method, values_from=tv) |>
    group_by(outcome) |>
    summarize(max_improvement = max(1 - model / weight),
              min_improvement = min(1 - model / weight))

filter(disp_party, level=="block", party=="dem") |>
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
ggsave(here("paper/figures/nc_overview.pdf"), plot=p, width=7.5, height=2.5)

# Disparity plots -----

lbl_race = function(race) {
    str_c(races[race], " (", percent(as.numeric(p_r[race]), 0.1), ")")
}

geos = c(county="County", zip="ZIP code", tract="Census tract", block="Census block")
geos_short = c(county="County", zip="ZIP", tract="Tract", block="Block")
methods = c(model="Model", ols="OLS", thresh="Threshold", weight="Weighted")

make_disp_plot = function(d, x, y, title, xlab) {
    filter(d, level=="block") |>
    ggplot(aes({{ x }}, {{ y }},
               color=methods[method], shape=methods[method])) +
        # facet_wrap(~ factor(geos[level], levels=geos)) +
        geom_hline(yintercept=0.0, color="#00000077") +
        geom_point(size=2.4, position=position_dodge(width=0.5)) +
        scale_y_continuous("Error in disparity estimation",
                           labels=label_number(1, scale=100, suffix="pp")) +
        labs(x=xlab, title=title, color="Method", shape="Method") +
        theme_paper()
}

p1 = make_disp_plot(disp_party, str_to_upper(party), disp_wb - disp_wb_true,
                    "White-Black disparity error (party)", "Party")
p2 = make_disp_plot(disp_party, str_to_upper(party), disp_wh - disp_wh_true,
                    "White-Hispanic disparity error (party)", "Party")
p3 = make_disp_plot(disp_turnout, fct_inorder(n_voted), disp_wb - disp_wb_true,
                    "White-Black disparity error (turnout)", "Number of elections voted")
p4 = make_disp_plot(disp_turnout, fct_inorder(n_voted), disp_wh - disp_wh_true,
                    "White-Hispanic disparity error (turnout)", "Number of elections voted")

p = p1 + p2 + p3 + p4 + plot_layout(guides="collect")
ggsave(here("paper/figures/nc_disp.pdf"), plot=p, width=8, height=5)


# TV plots -----

## Party fit quality -----
p1 = filter(tv_party, race=="overall") %>%
ggplot(aes(x=factor(geos[level], levels=geos), y=tv,
           color=methods[method], shape=methods[method], group=methods[method])) +
    geom_textline(aes(label=methods[method]),
                  position=position_dodge(width=0.25),
                  linewidth=0.7, size=3.5, family="Times", hjust=0.08) +
    geom_point(size=2.4, position=position_dodge(width=0.25)) +
    scale_color_wa_d() +
    scale_x_discrete(expand=c(0.07, 0, 0.06, 0)) +
    scale_y_log10("Overall total variation distance") +
    labs(x="BISG geographic precision", title="Party identification") +
    guides(color="none", shape="none") +
    theme_paper() +
    theme(plot.margin=unit(c(0, 0.1, 0, 0), "cm"))

p1b = filter(tv_party, race!="overall") %>%
ggplot(aes(x=factor(geos_short[level], levels=geos_short), y=tv,
           color=methods[method], shape=methods[method], group=methods[method])) +
    facet_wrap(~ factor(lbl_race(race), levels=lbl_race(names(races)))) +
    geom_line(position=position_dodge(width=0.25), size=0.7) +
    geom_point(size=2.0, position=position_dodge(width=0.25)) +
    scale_color_wa_d() +
    scale_y_log10("Total variation distance") +
    labs(x="BISG geographic precision", title="Party") +
    guides(color="none", shape="none") +
    theme_paper() +
    theme(plot.margin=unit(c(0, 0, 0, 0), "cm"))


## Turnout fit quality -----
p2 = filter(tv_turnout, race=="overall") %>%
ggplot(aes(x=factor(geos[level], levels=geos), y=tv,
           color=methods[method], shape=methods[method], group=methods[method])) +
    geom_textline(aes(label=methods[method],
                      hjust=c(model=0.08, weight=0.08, thresh=0.52, ols=0.1)[method]),
                  position=position_dodge(width=0.25),
                  linewidth=0.7, size=3.5, family="Times") +
    geom_point(size=2.4, position=position_dodge(width=0.25)) +
    scale_color_wa_d() +
    scale_x_discrete(expand=c(0.07, 0, 0.06, 0)) +
    scale_y_log10("Overall total variation distance") +
    labs(x="BISG geographic precision", title="Turnout") +
    theme_paper() +
    theme(plot.margin=unit(c(0, 0.1, 0, 0), "cm"))

p2b = filter(tv_turnout, race!="overall") %>%
ggplot(aes(x=factor(geos_short[level], levels=geos_short), y=tv,
           color=methods[method], shape=methods[method], group=methods[method])) +
    facet_wrap(~ factor(lbl_race(race), levels=lbl_race(names(races)))) +
    geom_line(position=position_dodge(width=0.25), size=0.7) +
    geom_point(size=2.0, position=position_dodge(width=0.25)) +
    scale_color_wa_d() +
    scale_y_log10("Total variation distance") +
    labs(x="BISG geographic precision", title="Turnout") +
    theme_paper() +
    theme(plot.margin=unit(c(0, 0, 0, 0), "cm"))

p = p1 + p2
ggsave(here("paper/figures/nc_tv.pdf"), plot=p, width=8, height=3.75)

p = p1b + p2b + plot_layout(ncol=1, guides="collect")
ggsave(here("paper/figures/nc_tv_detailed.pdf"), plot=p, width=6.5, height=8)
