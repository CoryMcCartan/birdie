suppressMessages({
    library(tidyverse)
    devtools::load_all(".")
    library(here)
})

d = readRDS(here("data-raw/nc_voters_small.rds")) |>
    filter(!is.na(party))
p_r = with(d, prop.table(table(race)))
p_y = with(d, prop.table(table(party)))
r_probs = bisg(~ nm(last_name) + zip(zip), data=d, p_r=p_r)

ctrl = birdie.ctrl(abstol=1e-3, reltol=1e-3, max_iter=100)
fit0 = birdie(r_probs, party ~ proc_zip(zip), data=d, ctrl=ctrl)
fit = birdie(r_probs, party ~ (1 | proc_zip(zip)), data=d, family=cat_mixed(), ctrl=ctrl)
fit2 = birdie(r_probs, party ~ (1 | proc_zip(zip)), data=d, family=cat_mixed(), ctrl=ctrl)

plot(t(coef(fit, subgroup=TRUE)[1,1:2,]), cex=0.5, xlim=0:1, ylim=0:1)
plot(t(coef(fit2, subgroup=TRUE)[1,1:2,]), cex=0.5, xlim=0:1, ylim=0:1)
fit
