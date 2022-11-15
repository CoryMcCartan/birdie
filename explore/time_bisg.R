library(microbenchmark)
devtools::load_all()
data("pseudo_vf")

# make tables
{
p_r = prop.table(table(pseudo_vf$race))

G_vec = as.factor(coalesce(pseudo_vf$zip, "<none>"))
p_rgx = census_zip_table(G_vec, "zip", p_r) |>
    rbind(rlang::list2(zip="<none>",
                       white=p_r[1], black=p_r[2], hisp=p_r[3],
                       asian=p_r[4], aian=p_r[5], other=p_r[6]))
p_rgx = mutate(p_rgx, across(.data$white:.data$other, ~ coalesce(., p_r[cur_column()])))
match_idx = match(levels(G_vec), p_rgx$zip)
match_idx[is.na(match_idx)] = match("<none>", p_rgx$zip)
p_rgx = p_rgx[match_idx, ]
p_gx = prop.table(table(G_vec))
p_gxr = as.matrix(select(p_rgx, .data$white:.data$other))
for (i in seq_along(p_r)) {
    p_gxr[, i] = p_gxr[, i] * p_gx
    p_gxr[, i] = p_gxr[, i] / sum(p_gxr[, i])
}

S_vec = as.character(pseudo_vf$last_name)
p_sr = census_surname_table(S_vec, "last_name", p_r, flip=TRUE)
S_vec[!S_vec %in% p_sr[[1]]] = "<generic>"
S_vec = factor(S_vec, levels=p_sr[[1]])
p_sr = as.matrix(p_sr[, -1])
}

S = rep(as.integer(S_vec), 100)
GX = rep(as.integer(G_vec), 100)


microbenchmark(
    est_bisg(S, GX, p_sr, p_gxr, p_r),
    mean(S),
    colMeans(p_sr[S, ])
    # est_bisg2(S, GX, p_sr, p_gxr, p_r)
)
