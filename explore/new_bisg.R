data(pseudo_vf)
p_r = prop.table(table(pseudo_vf$race))

formula = ~ nm(last_name) + zip(zip)

bisg(race ~ nm(last_name) + zip(zip), data=pseudo_vf)

bisg(race ~ nm(last_name), data=pseudo_vf)


# comparison

rpr1 = bisg(~ nm(last_name) + zip(zip), data=pseudo_vf, p_r=p_r)
rpr2 = predict_race_sgz(last_name, zip, data=pseudo_vf, p_r=p_r,
                        est_r_gz=FALSE, iterate=0)

hist(as.matrix(rpr1 - rpr2))

microbenchmark::microbenchmark(
    bisg(~ nm(last_name) + zip(zip), data=pseudo_vf),
    predict_race_sgz(last_name, zip, data=pseudo_vf, p_r=p_r,
                     est_r_gz=FALSE, iterate=0),
    times=10
)
