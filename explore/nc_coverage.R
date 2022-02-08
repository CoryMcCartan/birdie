library(tidyverse)
library(ggrepel)
library(scales)
library(wacolors)
suppressMessages(library(here))
devtools::load_all(here("."))

if (file.exists(voterfile <- here("data/nc_voters.rds"))) {
    voters = readRDS(voterfile)
} else {
    counties = tigris::fips_codes$county[tigris::fips_codes$state == "NC"] %>%
        stringr::str_sub(end=-8)
    voters = do.call(rbind, lapply(counties, make_nc_df))
    voters = voters %>%
        select(last_name, party, race, zip, gender, age, birth_state, lic) %>%
        mutate(across(where(is.character), factor))
    saveRDS(voters, voterfile, compress="xz")
}

voters = slice_sample(voters, n=1024) %>%
    mutate(gender = as.factor(coalesce(gender, "U"))) %>%
    filter(!is.na(age))

p_xr_true = prop.table(table(voters$party, voters$race))
p_r_true = colSums(p_xr_true)

run_boot = function(i) {
    d_fit = slice_sample(voters, n=nrow(voters), replace=T)

    p_xr = prop.table(table(d_fit$party, d_fit$race))
    p_r = colSums(p_xr)

    capture.output(suppressMessages({
        fit_bisg = model_race(party, last_name, zip, data=d_fit, p_r=p_r,
                              regularize=T, alpha=3, tol_rhat=1.25, lr=0.15,
                              methods=c("bis", "bisg", "pyro"))
        xr_bisg = calc_joints(p_xr, d_fit, fit_bisg)

        fit_true = model_race(party, last_name, zip, data=d_fit, p_r=p_r,
                              alpha=3, tol_rhat=1.25, lr=0.15,
                              methods=c("bis", "bisg", "pyro"),
                              use_true_gz=TRUE)
        xr_true = calc_joints(p_xr, d_fit, fit_true)
    }))

    tibble(draw = i,
           race = rep(levels(voters$race), each=4),
           party = rep(levels(voters$party), 5),
           true = as.numeric(p_xr_true %*% diag(1/p_r_true)),
           true_boot = as.numeric(p_xr %*% diag(1/p_r)),
           post_bisg_med = as.numeric(xr_bisg$pyro %*% diag(1/p_r)),
           post_true_med = as.numeric(xr_true$pyro %*% diag(1/p_r)),
           ci_width_bisg = as.numeric((xr_bisg$pyro_high - xr_bisg$pyro_low) %*% diag(1/p_r)),
           ci_width_true = as.numeric((xr_true$pyro_high - xr_true$pyro_low) %*% diag(1/p_r)),
           covered_bisg = as.numeric(xr_bisg$pyro_low < p_xr_true & p_xr_true < xr_bisg$pyro_high),
           covered_true = as.numeric(xr_true$pyro_low < p_xr_true & p_xr_true < xr_true$pyro_high),
           covered_bisg_boot = as.numeric(xr_bisg$pyro_low < p_xr & p_xr < xr_bisg$pyro_high),
           covered_true_boot = as.numeric(xr_true$pyro_low < p_xr & p_xr < xr_true$pyro_high))
}

R = 100
res = map_dfr(cli::cli_progress_along(1:R), run_boot) %>%
    mutate(race = factor(race, levels=levels(voters$race)),
           party = factor(party, levels=levels(voters$party)))

d = res %>%
    mutate(err_bisg = post_bisg_med - true,
           err_true = post_true_med - true) %>%
    group_by(party, race) %>%
    summarize(lbl = str_glue("{race[1]} {party[1]}"),
              ci_width_bisg = mean(ci_width_bisg),
              ci_width_true = mean(ci_width_true),
              bias_bisg = mean(err_bisg),
              bias_true = mean(err_true),
              std_bisg = sd(err_bisg),
              std_true = sd(err_true),
              coverage_bisg = mean(covered_bisg),
              coverage_true = mean(covered_true))

ggplot(d, aes(bias_bisg, coverage_bisg, label=lbl)) +
    geom_vline(xintercept=0, col="#00000077") +
    geom_hline(yintercept=0.9, col="#002233cc", lty="dashed") +
    geom_point(aes(shape=party, color=race), size=5) +
    geom_text_repel(size=2.8, fontface="bold", point.padding=6) +
    lims(x=c(-0.3, 0.35)) +
    scale_y_continuous("Coverage", labels=percent) +
    scale_fill_wa_d("larch")
ggplot(d, aes(bias_true, coverage_true, label=lbl)) +
    geom_vline(xintercept=0, col="#00000077") +
    geom_hline(yintercept=0.9, col="#002233cc", lty="dashed") +
    geom_point(aes(shape=party, color=race), size=5) +
    geom_text_repel(size=2.8, fontface="bold", point.padding=6) +
    lims(x=c(-0.3, 0.35)) +
    scale_y_continuous("Coverage", labels=percent) +
    scale_fill_wa_d("larch")

ggplot(d, aes(bias_bisg, std_bisg, label=lbl)) +
    geom_vline(xintercept=0, col="#00000077") +
    geom_point(aes(shape=party, color=race), size=5) +
    geom_text_repel(size=2.8, fontface="bold", point.padding=6) +
    scale_y_continuous("Std. dev. in error of posterior median", trans="log10")

ggplot(d, aes(std_bisg, ci_width_bisg, label=lbl)) +
    geom_abline(slope=1, col="#00000077") +
    geom_point(aes(shape=party, color=race), size=5) +
    geom_text_repel(size=2.8, fontface="bold", point.padding=6) +
    scale_x_continuous("Std. dev. in error of posterior median", trans="log10", limits=c(0.002, 0.3)) +
    scale_y_continuous("Mean credible interval width", trans="log10", limits=c(0.007, 0.8))
ggplot(d, aes(std_true, ci_width_true, label=lbl)) +
    geom_abline(slope=1, col="#00000077") +
    geom_point(aes(shape=party, color=race), size=5) +
    geom_text_repel(size=2.8, fontface="bold", point.padding=6) +
    scale_x_continuous("Std. dev. in error of posterior median", trans="log10", limits=c(0.002, 0.3)) +
    scale_y_continuous("Mean credible interval width", trans="log10", limits=c(0.007, 0.8))
ggsave("~/Desktop/ss.png", width=9.5, height=4.5, dpi=600)

ggplot(d, aes(coverage_bisg, coverage_true, label=lbl)) +
    geom_hline(yintercept=0.9, col="#002233cc", lty="dashed") +
    geom_vline(xintercept=0.9, col="#002233cc", lty="dashed") +
    geom_point(aes(shape=party, color=race), size=5) +
    geom_text_repel(size=2.8, fontface="bold", point.padding=6) +
    scale_x_continuous("Coverage", labels=percent) +
    scale_y_continuous("Coverage", labels=percent, limits=c(0.35, 1))

ggplot(d, aes(bias_bisg, bias_true, label=lbl)) +
    geom_point(aes(shape=party, color=race), size=5) +
    geom_text_repel(size=2.8, fontface="bold", point.padding=6)


if (F) {
res %>%
    filter(!race %in% c("asian", "other")) %>%
    mutate(race = factor(race, levels=levels(voters$race)),
           party = factor(party, levels=levels(voters$party)),
           lbl = str_glue("{race} {party}"),
           err = post_med - true,
           coverage = post_med - ci_width/2 < true &
               post_med + ci_width/2 > true) %>%
ggplot(aes(race, err, fill=race)) +
    geom_hline(yintercept=0, col="#00000077") +
    facet_wrap(~ party, scales="free_y") +
    geom_boxplot() +
    scale_y_continuous("Posterior median - true", trans="pseudo_log")

}


run_boot = function(i) {
    d_fit = slice_sample(voters, n=nrow(voters), replace=T)
    p_r = prop.table(table(d_fit$race))

    fit = model_race(party, last_name, zip, data=d_fit, p_r=p_r,
                     regularize=T, alpha=3, tol_rhat=1.25, lr=0.15, iter=500,
                     methods=c("bisg", "pyro"))

    fit$pyro$p_xr
}

fit_boot = model_race(party, last_name, zip, data=d_fit, p_r=p_r,
                      regularize=T, alpha=3,
                      methods=c("bis", "bisg"))
res_boot = map(cli::cli_progress_along(1:10), run_boot)
fit_boot$pyro = list(p_xr = do.call(abind::abind, c(res_boot, list(along=1))))

xr = calc_joints(p_xr, voters, fit_boot)
