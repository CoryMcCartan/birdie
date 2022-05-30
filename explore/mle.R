suppressMessages({
    library(here)
    devtools::load_all(here())
})

data(pseudo_vf)

p_xr = prop.table(table(pseudo_vf$turnout, pseudo_vf$race))
p_r = prop.table(table(pseudo_vf$race))

r_probs = predict_race_sgz(last_name, zip, data=pseudo_vf, p_r=p_r)

fit_vi = model_race(r_probs, turnout, G=zip, data=pseudo_vf, method="vi")
fit_mle = model_race(r_probs, turnout, G=zip, data=pseudo_vf, method="mle")

xr = list(
    true = p_xr,
    vi = calc_joint_model(fit_vi, p_r=p_r),
    mle = calc_joint_model(fit_mle, p_r=p_r)
)

eval_joints(xr$true, "tv",
            vi = xr$vi,
            mle = xr$mle)
