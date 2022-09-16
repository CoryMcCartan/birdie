library(tidybayes)
library(posterior)

voters = read_rds(here("data-raw/nc_voters.rds")) |>
    filter(county_name == "Dare") |>
    slice_sample(n=2000, replace=FALSE)

r_probs = predict_race_sgz(last_name, zip, data=voters)

N_fit = 100
N_draws = 1000
fit_svi = model_race(head(r_probs, N_fit), party, zip,
                     data=head(voters, N_fit), config=list(lr=0.3, n_draws=N_draws))

fit_hmc = model_race_hmc(head(r_probs, N_fit), party, zip,
                         data=head(voters, N_fit), config=list(warmup=500, n_draws=N_draws))


d_draws = tibble(.iter = rep(seq_len(N_draws), 4*6),
       party = rep(rep(fit_hmc$x_lev, each=N_draws), 6),
       race = rep(fit_hmc$r_lev, each=N_draws*4),
       draw_svi = as.numeric(fit_svi$draws$global),
       draw_nuts = as.numeric(fit_hmc$draws$global)) |>
    pivot_longer(draw_svi:draw_nuts, names_prefix="draw_", names_to="alg")

ggplot(d_draws, aes(toupper(party), value, fill=toupper(alg))) +
    facet_wrap(~ races[race]) +
    geom_boxplot(outlier.size=0.2) +
    labs(x="Party", y="Estimate", fill="Algorithm") +
    scale_fill_wa_d("baker", which=c(2, 12)) +
    theme_paper()

ggsave(here("paper/figures/valid_nuts_compare.pdf"), width=7, height=4.5)
