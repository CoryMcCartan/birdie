library(tidyverse)
library(ggrepel)
library(scales)
library(wacolors)
suppressMessages(library(here))
devtools::load_all(here("."))

suppressMessages(library(here))
devtools::load_all(here("."))

if (file.exists(voterfile <- here("data/nc_voters.rds"))) {
    voters = readRDS(voterfile)
} else {
    voters = make_nc_statewide(voterfile)
}

voters = slice_sample(voters, n=1024) %>%
    mutate(gender = as.factor(coalesce(gender, "U"))) %>%
    filter(!is.na(age))


# tables
p_xr_true = prop.table(table(voters$party, voters$race))
p_r_true = colSums(p_xr_true)


run_boot = function(i) {
    d_fit = slice_sample(voters, n=nrow(voters), replace=T)

    p_xr = prop.table(table(d_fit$party, d_fit$race))
    p_r = colSums(p_xr)


    capture.output(suppressMessages({
        r_probs = predict_race_sgz(last_name, zip, c(gender, age), data=d_fit,
                                   p_r=p_r, iterate=0)

        fit0 = model_race(r_probs$pr, party, zip, c(gender, age), data=d_fit, sgz=r_probs,
                          reload_py=F, config=list(n_mi=0, lr=0.2), silent=T)
        fit1 = model_race(r_probs$pr, party, zip, c(gender, age), data=d_fit, sgz=r_probs,
                          reload_py=F, config=list(n_mi=5, lr=0.2), silent=T)

        xr = list(
            true = p_xr,
            bisgz = calc_joint_bisgz(r_probs$pr, d_fit$party),
            mod = calc_joint_model(fit0$p_xr, 0.5, p_r),
            mod_low = calc_joint_model(fit0$p_xr, 0.05, p_r),
            mod_high = calc_joint_model(fit0$p_xr, 0.95, p_r),
            mi = calc_joint_model(fit1$p_xr, 0.5, p_r),
            mi_low = calc_joint_model(fit1$p_xr, 0.05, p_r),
            mi_high = calc_joint_model(fit1$p_xr, 0.95, p_r)
        )
    }))

    tibble(draw = i,
           race = rep(levels(voters$race), each=4),
           party = rep(levels(voters$party), 5),
           true = as.numeric(p_xr %*% diag(1/p_r)),
           post_mod_med = as.numeric(xr$mod %*% diag(1/p_r)),
           post_mi_med = as.numeric(xr$mi %*% diag(1/p_r)),
           ci_width_mod = as.numeric((xr$mod_high - xr$mod_low) %*% diag(1/p_r)),
           ci_width_mi = as.numeric((xr$mi_high - xr$mi_low) %*% diag(1/p_r)),
           covered_mod = as.numeric(xr$mod_low < p_xr & p_xr < xr$mod_high),
           covered_mi = as.numeric(xr$mi_low < p_xr & p_xr < xr$mi_high))
}

R = 100
res = map_dfr(cli::cli_progress_along(1:R), run_boot) %>%
    mutate(race = factor(race, levels=levels(voters$race)),
           party = factor(party, levels=levels(voters$party)))

d = res %>%
    mutate(err_mod = post_mod_med - true,
           err_mi = post_mi_med - true,
           cover_diff = covered_mi - covered_mod) %>%
    group_by(party, race) %>%
    summarize(lbl = str_glue("{race[1]} {party[1]}"),
              ci_width_mod = mean(ci_width_mod),
              ci_width_mi = mean(ci_width_mi),
              bias_mod = mean(err_mod),
              bias_mi = mean(err_mi),
              std_mod = sd(err_mod),
              std_mi = sd(err_mi),
              coverage_mod = mean(covered_mod),
              coverage_mi = mean(covered_mi),
              coverage_diff = mean(cover_diff))

ggplot(d, aes(bias_mod, coverage_diff, label=lbl)) +
    geom_vline(xintercept=0, col="#00000077") +
    #geom_hline(yintercept=0.9, col="#002233cc", lty="dashed") +
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

ci_mult = diff(qnorm(c(0.05, 0.95)))
ggplot(d, aes(ci_mult*std_mi, ci_width_mi, label=lbl)) +
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
