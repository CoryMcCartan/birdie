devtools::load_all(".")
data(pseudo_vf)
pseudo_vf = pseudo_vf |>
    mutate(y = 1*(turnout == "yes"))

form = turnout ~ 0 + white + black + hisp + asian + aian + other

races = levels(pseudo_vf$race)
d_race_pr = as_tibble(as.data.frame(diag(6)))
colnames(d_race_pr) = races

m_race = model.matrix(~ 0 + race, data=pseudo_vf)
colnames(m_race) = races
m_true = glm(form, bind_cols(pseudo_vf, m_race), family=binomial())

get_pred <- function(m) {
    out = predict(m, newdata=d_race_pr, type="response")
    names(out) = races
    out
}

px = mean(pseudo_vf$turnout == "yes")
xr_act = get_pred(m_true)

yes = which(pseudo_vf$turnout == "yes")

r_probs = predict_race_sgz(last_name, zip, data=pseudo_vf)
# r_probs2 = predict_race_sgz_me(last_name, zip, data=pseudo_vf)
# r_probs = r_probs2
m_race = as.matrix(r_probs)
colnames(m_race) = races
m_race0 = m_race

xr_ols = coef(lm(turnout=="yes" ~ 0 + white + black + hisp + asian + aian + other,
                 data=bind_cols(pseudo_vf, m_race0)))

# xr_est = xr_ols
xr_est = colSums(m_race0[yes, ])/colSums(m_race0)
# xr_est = rbind(px * colMeans(m_race[yes, ]), (1 - px) * colMeans(m_race[-yes, ]))
# xr_est = (xr_est %*% diag(1/colSums(xr_est)))[1, ]
ests = xr_est

for (i in 1:25) {
    # E step
    m_race[yes, ] = m_race0[yes, ] %*% diag(xr_est)
    m_race[-yes, ] = m_race0[-yes, ] %*% diag(1 - xr_est)
    m_race = m_race / rowSums(m_race)
    # plot(m_race0[, 1], m_race[, 1], cex=0.01)

    # M step
    # xr_est = rbind(px * colMeans(m_race[yes, ]), (1 - px) * colMeans(m_race[-yes, ]))
    xr_est = colSums(m_race[yes, ])/colSums(m_race)
    # for (j in 1:6) {
    #     x = lme4::lmer(I(turnout == "yes") ~ (1 | zip), data=pseudo_vf,
    #                     weights=m_race[, j]) |>
    #         suppressMessages() |>
    #         suppressWarnings()
    #     xr_est[j] = x@beta#plogis(x@beta)
    # }
    ests = rbind(ests, xr_est)
    # xr_est = (xr_est %*% diag(1/colSums(xr_est)))[1, ]
    # m = glm(form, bind_cols(pseudo_vf, m_race), family=binomial())
    # xr_est = get_pred(m)*0.1 + 0.9*xr_est
    # xr_est = get_pred(m)
}
matplot(ests, type='b', cex=0.5)

round(xr_est, 3)
sum(abs(xr_est[1:5] - xr_act[1:5]))
sum(abs(xr_ols[1:5] - xr_act[1:5]))

pseudo_vf$dummy = sample(c("A", "B"), nrow(pseudo_vf), replace=TRUE)
fit = model_race(r_probs, turnout, dummy, data=pseudo_vf)
xr_mod = colMeans(fit$draws$global)[2,]
xr_mod
