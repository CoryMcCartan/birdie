library(matrixStats)
library(collapse)

d_surname <- read_rds(here("data-out/surn_cov.rds")) |>
    mutate(nm_grp = fct_na_value_to_level(fct_collapse(
        group,
        anglobl = c("anglo", "black"),
        eur_wave1 = c("german", "nordic", "irel"),
        eur_wave2 = c("eeur", "italy", "jwsrus", "agean"),
        easia = c("china", "japkor"),
        seasiapac = c("seasia", "asiapac"),
        hisp = c("mexico", "prq", "latin"),
        other = c("weur", "mena", "native", "afrocar")
    ), "other"))
print(table(d_surname$group))
print(table(d_surname$nm_grp))

# join to data and do sensitivity analysis ########

d <- d |>
    left_join(d_surname, by="last_name") |>
    mutate(nm_grp = fct_na_value_to_level(nm_grp, "other"))

m_sens = birdie(r_probs$county, party ~ nm_grp + GEOID_county, data=d, ctrl=ctrl) |>
    suppressWarnings()

p_xr = prop.table(table(d$party, d$race), 2)
diffs = coef(m_sens) - coef(fits_party$county$sat)
write_rds(diffs, here("paper/data/nc_sens_diff.rds"))

colSums(abs(coef(m_sens) - p_xr)) %*% p_r * 0.5
colSums(abs(coef(fits_party$county$sat) - p_xr)) %*% p_r * 0.5

res = residuals(fits_party$county$sat)

calc_cov = function(grp, n=10000) {
    m_cov = cov(cbind(d$nm_grp == grp, res))
    qs = rowQuantiles(
        rWishart(n, nrow(d), m_cov)[1,-1,] / nrow(d),
        probs=c(0.05, 0.95)
    )
    scl = sqrt(m_cov[1, 1] * diag(m_cov)[-1])
    qTBL(data.frame(group=toupper(grp), party=toupper(colnames(res)),
                    cov=m_cov[-1, 1]/scl, cov_q05=qs[, 1]/scl, cov_q95=qs[, 2]/scl))
}

d_cov = levels(d$nm_grp) |>
    map_dfr(calc_cov) |>
    mutate(sign = sign(cov) == sign(cov_q05) & sign(cov) == sign(cov_q95),
           group = fct_relevel(group,
                               "ANGLOBL", "EUR_WAVE1", "EUR_WAVE2",
                               "HISP", "CUBA", "EASIA", "SEASIAPAC", "SASIA"))

ggplot(d_cov, aes(group, fct_inorder(party), fill=cov)) +
    geom_tile() +
    geom_text(aes(label=if_else(sign, "*", "")), size=4, family="Times") +
    coord_cartesian(expand=F) +
    scale_fill_wa_c("vantage", midpoint=0, labels=percent, reverse=TRUE) +
    labs(y="Party", x="Surname group", fill="Residual\ncorrelation") +
    theme_paper()

ggsave(here("paper/figures/nc_sens_grp.pdf"), width=7.5, height=3.5)


bind_rows(
    sens=tidy(m_sens),
    orig=tidy(fits_party$county$sat),
    .id="fit"
) |>
    pivot_wider(names_from=fit, values_from=estimate) |>
ggplot(aes(orig, sens, color=toupper(party))) +
    # geom_hline(yintercept=0, lty="dashed") +
    geom_abline(slope=1, lty="dashed", color="#444444") +
    facet_wrap(~ race) +
    geom_point(size=3) +
    scale_color_manual(values=c(LIB="#E4C22B", REP="#B25D4C", IND="#aaccaf", DEM="#3D77BB")) +
    theme_paper()

# save names for appendix ########

d_surname |>
    group_by(nm_grp) |>
    slice_head(n=50) |>
    select(-group) |>
    ungroup() |>
write_rds("paper/data/nm_grps.rds", compress="xz")

X = Matrix::sparse.model.matrix(~ last_name, data=d)
m_l = glmnet::glmnet(X, rowMeans(abs(res)))
which(coef(m_l, s=m_l$lambda[10])[-1, 1] != 0)
ff = predict(m_l, s=m_l$lambda[10], newx=X)[,1]
cor(ff, res)
