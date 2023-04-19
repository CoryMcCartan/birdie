source("replication/00_setup.R")

# Setup --------
voters = read_rds(here("data-raw/nc_voters.rds"))
set.seed(5118)
d = voters |>
    mutate(GEOID_county = as.character(county),
           GEOID_tract = if_else(is.na(tract), GEOID_county, str_c(county, tract)),
           GEOID_block = if_else(is.na(block), GEOID_county, str_c(county, tract, block)),
           GEOID_zip = if_else(is.na(zip), str_c("cty", county), as.character(zip)),
           party=coalesce(party, "ind")) |>
    slice_sample(n=1e6)
print(as.character(head(d$last_name))) # ensure seed is working


# Do BISG --------
p_r = prop.table(table(d$race))
r_probs = bisg(~ nm(last_name), data=d, p_r=p_r)

geo_levels = c("county", "zip", "tract")
ctrl = birdie.ctrl(abstol=1e-5)
fits_party <- map(geo_levels, function(level) {
    cat(level, "\n")

    d_mm = rename(d, GEOID=str_c("GEOID_", level))

    suppressWarnings(list(
        sat = birdie(r_probs, party ~ GEOID, data=d_mm, ctrl=ctrl)
    ))
}) |>
    set_names(geo_levels)

# run functions at top of 02_nc_fig.R

r_probs = list(county=r_probs, zip=r_probs, tract=r_probs)
ests_party = make_est_d(fits_party, party, d)

disp_party = calc_disp(ests_party)

tv_party = eval_fit_tv(ests_party)


p1 = make_disp_plot(disp_party, str_to_upper(party), disp_wb - disp_wb_true,
                    "White-Black disparity error", "Party", lev="tract")
p2 = make_disp_plot(disp_party, str_to_upper(party), disp_wh - disp_wh_true,
                    "White-Hispanic disparity error", "Party", lev="tract")

p = p1 + p2 + plot_layout(guides="collect")



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
                       limits=c(0.0, NA), expand=expansion(c(0, 0.05))) +
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
    scale_y_continuous("Total variation distance", limits=c(0, NA),
                       expand=expansion(c(0, 0.03))) +
    labs(x="BISG geographic precision", #title="Party",
         color="Method", shape="Method") +
    scale_color_manual(values=methods_col, labels=methods) +
    scale_shape_manual(values=methods_shp, labels=methods) +
    guides(shape="none", color="none") +
    theme_paper() +
    theme(plot.margin=unit(c(0, 0, 0, 0), "cm"))


p = p1 + p1b
ggsave(here("~/Desktop/nc_tv.pdf"), plot=p, width=8, height=4)
