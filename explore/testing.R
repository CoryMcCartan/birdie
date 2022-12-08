suppressMessages({
    library(tidyverse)
    devtools::load_all(".")
    library(cluster)
    library(here)
})

d = readRDS(here("data-raw/nc_voters_small.rds")) |>
    filter(!is.na(party)) |>
    mutate(ln100 = fct_lump_n(last_name, 100))

p_r = with(d, prop.table(table(race)))
p_y = with(d, prop.table(table(party)))

d_rs = census_surname_table(d$last_name, "last_name")
p_rs = as.matrix(d_rs[-1])
cl = clara(p_rs, k=200, samples=100, cluster.only=TRUE)
# d_rs$last_name[cl == 10]

idx = match(d$last_name, d_rs$last_name)
ok = which(!is.na(idx))
X = model.matrix(~ 0 + as.factor(cl[idx]))

rr = resid(lm(n_voted == "5" ~ race, data=d[ok, ]))
m = lm.fit(X, rr)
cor(m$fitted.values, rr)^2

rr = lm.fit(X, d$n_voted[ok] == "5")$residuals
m = lm(rr ~ 0 + race, data=d[ok, ])
summary(m)$r.squared

# grp = split(ok, fct_cross(d$race[ok], d$county[ok]))
grp = split(ok, cl[idx][ok])
res = map_dbl(grp, possibly(function(ii) {
    # chisq.test(cl[idx][ii], d$party[ii], simulate.p.value=TRUE)$p.value
    # chisq.test(d$last_name[ii], d$party[ii], simulate.p.value=TRUE, B=200)$p.value
    chisq.test(d$race[ii], d$party[ii], simulate.p.value=TRUE, B=1000)$p.value
}, NA))

