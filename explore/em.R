devtools::load_all(".")
data(pseudo_vf)

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

m_race = as.matrix(predict_race_sgz(last_name, zip, data=pseudo_vf))
colnames(m_race) = races
m_race0 = m_race

xr_ols = coef(lm(turnout=="yes" ~ 0 + white + black + hisp + asian + aian + other,
                 data=bind_cols(pseudo_vf, m_race0)))

xr_est = xr_ols
# xr_est = rbind(px * colMeans(m_race[yes, ]), (1 - px) * colMeans(m_race[-yes, ]))
# xr_est = (xr_est %*% diag(1/colSums(xr_est)))[1, ]
for (i in 1:5) {
    # E step
    m_race[yes, ] = m_race0[yes, ] %*% diag(xr_est)
    m_race[-yes, ] = m_race0[-yes, ] %*% diag(1 - xr_est)
    m_race = m_race / rowSums(m_race)
    plot(m_race0[, 1], m_race[, 1], cex=0.01)

    # M step
    # xr_est = rbind(px * colMeans(m_race[yes, ]), (1 - px) * colMeans(m_race[-yes, ]))
    # xr_est = (xr_est %*% diag(1/colSums(xr_est)))[1, ]
    m = glm(form, bind_cols(pseudo_vf, m_race), family=binomial())
    xr_est = get_pred(m)*0.1 + 0.9*xr_est
    # xr_est = get_pred(m)
}

round(xr_est, 3)
sum(abs(xr_est - xr_act))
sum(abs(xr_ols - xr_act))
