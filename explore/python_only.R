data(pseudo_vf)

r_probs = predict_race_sgz(last_name, zip, data=pseudo_vf)

fit_objs = model_race_prep(r_probs, turnout, zip, data=pseudo_vf)

library(arrow)

write_parquet(fit_objs$X, "data-raw/python_only/X.parquet")
write_parquet(fit_objs$GZ_mat, "data-raw/python_only/GZ_mat.parquet")
write_parquet(fit_objs$GZ_var, "data-raw/python_only/GZ_var.parquet")
write_parquet(fit_objs$r_probs, "data-raw/python_only/r_probs.parquet")
write_parquet(fit_objs$preds, "data-raw/python_only/preds.parquet")

fit = reticulate::py_load_object("data-raw/python_only/fit.pkl")
class(fit) = "fit_raceproxy"
fit$N = nrow(r_probs)
fit$vars = "zip"
fit$x_lev = levels(pseudo_vf$turnout)
fit$r_lev = colnames(r_probs)

fit
