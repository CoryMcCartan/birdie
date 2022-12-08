# run replication through 01_nc_fig.R

party = "dem"

m_rxs = as.matrix(r_probs$tract)*0.999 + 0.001
m_rx = as.matrix(attr(r_probs$tract, "p_rgx")[attr(r_probs$tract, "gx"), ])*0.999 + 0.001

e_a = colMeans(1 / m_rxs)
e_as = colMeans(1 / m_rx)

c_r = sqrt((e_a - e_as) / e_as)

d_yrx = tidy(fits_party$tract$mmm, subgroup=TRUE) |>
    filter(.data$party == .env$party) |>
    select(race, GEOID, estimate) |>
    pivot_wider(names_from=race, values_from=estimate)
m_yrx = left_join(d, d_yrx, by=c("GEOID_tract"="GEOID")) |>
    select(white:other) |>
    as.matrix()
p_yrxs = rowSums(m_yrx * m_rxs)
p_yrx = rowSums(m_yrx * m_rx)

e_res = mean(((d$party == party) - p_yrx)^2)

# comparison r^2
summary(lm(party == "dem" ~ GEOID_tract, data=d))$r.squared # Y | R, G
cor(p_yrx, d$party == party)^2 # Y | R, G
cor(p_yrxs, d$party == party)^2 # Y | R, G, S

s = sqrt(e_res * e_as)

r2 = seq(0, 1, 0.01)^2

d_bounds = map_dfr(r2, function(rr) {
    est = coef(fits_party$tract$mmm)[party, ]
    lb = est * (1 - c_r * s * sqrt(rr))
    ub = est * (1 + c_r * s * sqrt(rr))
    tibble(r2 = rr,
           race = names(lb),
           lb = pmax(lb, 0),
           est = est,
           ub = pmin(ub, 1))
})

ggplot(d_bounds, aes(r2, est, ymin=lb, ymax=ub, fill=race, color=race)) +
    geom_ribbon(alpha=0.2, color=NA) +
    geom_line() +
    coord_cartesian(xlim=c(0, 0.1), expand=FALSE)

# disp bounds
e_a2 = e_a["pr_white"] + e_a["pr_black"]
e_as2 = e_as["white"] + e_as["black"]
c_r2 = sqrt((e_a2 - e_as2) / e_as2)
s2 = sqrt(e_res * e_as2)

d_bounds = map_dfr(r2, function(rr) {
    est = diff(coef(fits_party$tract$mmm)[party, c("white", "black")])
    lb = est * (1 - c_r2 * s2 * sqrt(rr))
    ub = est * (1 + c_r2 * s2 * sqrt(rr))
    tibble(r2 = rr,
           lb = lb,
           est = est,
           ub = ub)
})

ggplot(d_bounds, aes(r2, est, ymin=lb, ymax=ub)) +
    geom_ribbon(alpha=0.25, color=NA) +
    geom_line() +
    coord_cartesian(ylim=c(-1, 1), xlim=c(0, 0.2), expand=FALSE)
