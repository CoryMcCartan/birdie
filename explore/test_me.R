library(raceproxy)
library(dplyr)

if (file.exists(voterfile <- here("data/nc_voters.rds"))) {
    voters = readRDS(voterfile)
} else {
    voters = make_nc_statewide(voterfile)
}

d_fit = slice_sample(voters, n=200e3) %>%
    mutate(gender = as.factor(coalesce(gender, "U"))) %>%
    filter(!is.na(age))

p_xr = prop.table(table(d_fit$party, d_fit$race))
p_x = rowSums(p_xr)
p_r = colSums(p_xr)

r_probs = predict_race_sgz(last_name, zip, data=d_fit, p_r=p_r, regularize=F, iterate=0)
r_probs_me = predict_race_sgz_me(last_name, zip, data=d_fit, p_r=p_r)

plot(as.matrix(r_probs)[,1], as.matrix(r_probs_me)[,1], cex=0.05)
plot(as.matrix(r_probs)[,2], as.matrix(r_probs_me)[,2], cex=0.05)
plot(as.matrix(r_probs)[,3], as.matrix(r_probs_me)[,3], cex=0.05)
plot(as.matrix(r_probs)[,4], as.matrix(r_probs_me)[,4], cex=0.05)
plot(as.matrix(r_probs)[,5], as.matrix(r_probs_me)[,5], cex=0.05)
plot(as.matrix(r_probs)[,6], as.matrix(r_probs_me)[,6], cex=0.05)

fit_no = model_race(r_probs, party, zip, data=d_fit)
fit_me = model_race(r_probs_me, party, zip, data=d_fit)

xr = list(
    true = p_xr,
    bisg_no = calc_joint_bisgz(r_probs, d_fit$party),
    bisg_me = calc_joint_bisgz(r_probs_me, d_fit$party),
    model_no = calc_joint_model(fit_no$p_xr, 0.5, p_r),
    model_me = calc_joint_model(fit_me$p_xr, 0.5, p_r)
)

eval_joints(xr$true, "tv",
            bisg_no = xr$bisg_no,
            bisg_me = xr$bisg_me,
            model_no = xr$model_no,
            model_me = xr$model_me)


races = c("White", "Black", "Hispanic", "Asian", "Native", "Other")
d_cond = xr[c("true", "model_no", "model_me")] %>%
    map_dfr(function(m) {
        tibble(pr = as.numeric(m %*% diag(1/p_r)),
               party = rep(levels(d_fit$party), 6),
               race = rep(races, each=4))
    }, .id="version") %>%
    mutate(version = factor(c(true="True", model_no="Orig.", model_me="ME")[version],
                            levels=c("ME", "True", "Orig.")),
           race = factor(race, levels=races),
           party = factor(str_to_upper(party),
                          levels=str_to_upper(levels(d_fit$party))))

ggplot(d_cond, aes(version, pr, fill=party)) +
    facet_grid(~ race) +
    geom_col(position=position_fill()) +
    scale_y_continuous("Proportion", labels=scales::percent) +
    scale_fill_manual(values=c(DEM="#3D77BB", IND="#aaccaf", REP="#B25D4C", LIB="#E4C22B")) +
    coord_cartesian(expand=F) +
    labs(x=NULL, fill="Party") +
    theme_minimal(base_family="Overpass", base_size=12) +
    theme(panel.spacing.x = unit(1, "lines"),
          strip.text=element_text(face="bold"))
ggsave("~/Desktop/out.png", width=11, height=5.1, dpi=400)
