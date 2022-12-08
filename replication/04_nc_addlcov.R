# Additional covariate --------

d$party_20 = fct_cross(d$party, d$voted_20)
ctrl = birdie.ctrl(abstol=1e-5)

suppressWarnings({
    fit_pt = birdie(r_probs$block, party_20 ~ (1 | GEOID_tract), data=d, ctrl=ctrl)
    fit_t = birdie(r_probs$block, voted_20 ~ (1 | GEOID_tract), data=d, ctrl=ctrl)
    fit_p = birdie(fitted(fit_t), party ~ voted_20 + (1 | GEOID_tract), data=d, ctrl=ctrl)
})

true_pt = tidy_true(d$race, d$party_20, "party_20") |>
    separate(party_20, c("party", "voted_20"), sep=":") |>
    group_by(voted_20, race) |>
    mutate(est_sum = sum(estimate)) |>
    ungroup() |>
    mutate(est_true = estimate / est_sum) |>
    select(-est_sum, -estimate)

est_pt_joint = tidy(fit_pt) |>
    separate(party_20, c("party", "voted_20"), sep=":") |>
    group_by(voted_20, race) |>
    mutate(est_sum = sum(estimate)) |>
    ungroup() |>
    mutate(est = estimate / est_sum,
           method="joint") |>
    select(-est_sum, -estimate)

est_pt_two = fitted(fit_p) |>
    est_weighted(party ~ voted_20, data=d) |> # marginalize out geo from `fit_t`
    tidy(subgroup=TRUE) |>
    rename(est = estimate) |>
    mutate(method="two", .before=everything())

est_pt_wtd = est_weighted(r_probs$block, party ~ voted_20, data=d) |>
    tidy(subgroup=TRUE) |>
    rename(est = estimate) |>
    mutate(method="weight", .before=everything())

idx_yes = which(d$voted_20 == "yes")
est_pt_thresh = bind_rows(
    yes = tidy_thresh(r_probs$block[idx_yes, ], d$party[idx_yes], nm="party"),
    no = tidy_thresh(r_probs$block[-idx_yes, ], d$party[-idx_yes], nm="party"),
    .id="voted_20"
) |>
    rename(est = estimate) |>
    mutate(method="thresh", .before=everything())

ests = bind_rows(est_pt_two, est_pt_joint, est_pt_wtd, est_pt_thresh) |>
    mutate(method = fct_inorder(method)) |>
    left_join(true_pt, by=c("voted_20", "party", "race")) |>
    mutate(level = "zip",
           vote_party = str_c(toupper(party), "/", str_to_title(voted_20)),
           .before=everything())

tv_ests = eval_fit_tv(ests)

l_tv = tv_ests |>
    filter(race == "overall") |>
    select(method, tv) |>
    pivot_wider(names_from=method, values_from=tv) |>
    as.list()

write_rds(l_tv, here("paper/data/nc_addlcov.rds"))

## plot -----

lbl_turn <- function(x) {
    c(yes="Voted in 2020", no="Did not vote in 2020")[x]
}

methods = c(two="BIRDiE (two-step)", joint="BIRDiE (joint)",
            thresh="Threshold", weight="Weighting")
methods_col = c(two=PAL_R[3], joint=PAL_R[3], thresh=PAL_R[4], weight=PAL_R[1])
methods_shp = c(two=16, joint=15, thresh=4, weight=1)

ggplot(ests, aes(str_to_upper(party), est - est_true, color=method, shape=method)) +
    facet_grid(fct_inorder(lbl_turn(voted_20)) ~ fct_inorder(lbl_race(race))) +
    geom_hline(yintercept=0, lty="dashed") +
    geom_point(size=2.0, position=position_dodge(width=0.7)) +
    scale_color_manual(values=methods_col, labels=methods) +
    scale_shape_manual(values=methods_shp, labels=methods) +
    scale_y_continuous("Error in probability estimation",
                       labels=label_number(1, scale=100, suffix="pp")) +
    guides(color=guide_legend(override.aes=aes(size=4))) +
    labs(x="Party", color="Method", shape="Method") +
    theme_paper() +
    theme(legend.position="bottom",
          legend.margin=margin())

ggsave(here("paper/figures/nc_addlcov_error.pdf"), width=8, height=4.5)



