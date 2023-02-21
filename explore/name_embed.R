library(dplyr)
library(purrr)
library(stringr)

d_names = readRDS(system.file("extdata", "names_2010_counts.rds",
                              package="birdie", mustWork=TRUE)) |>
    mutate(soundex = stringdist::phonetic(last_name, useBytes=TRUE))

cosine <- function(x, y) {
    sum(x * y) / sqrt(sum(x^2) + sum(y^2))
}
ngram_id <- function(str, n=2) {
    # if (novowel) str = str_remove_all(str, "(?<=.)[AEIOUY]")
    letters = as.integer(charToRaw(str)) - 65L
    seqn = seq_len(n) - 1
    windows = lapply(seq_len(length(letters) - n + 1), \(x) x + seqn)
    # vapply(windows, \(w) as.integer(sum(letters[w] * 26^seqn)), 0L)
    vapply(windows, \(w) rawToChar(as.raw(letters[w] + 65L)), "")
}

name_letters = d_names$last_name |>
    str_remove_all("(?<=.)[AEIOUY]") |>
    lapply(ngram_id, n=2)
ngram_tbl = sort(unique(unlist(name_letters)))
ngrams = name_letters |>
    map(\(x) tabulate(fastmatch::fmatch(x, ngram_tbl), length(ngram_tbl))) |>
    # map(\(x) x / (sum(x) + 1)) |>
    do.call(rbind, args=_)

sv <- as.matrix(dist(ngrams[1:1000, ])) |> svd(nu=5, nv=5)
pca <- prcomp(ngrams[1:10000, ], rank=10, center=T, scale.=F)
plot(pca$rotation[1:100, 2:3], type='n', xlim=c(-0.02, 0.02), ylim=c(-0.01, 0.01))
text(pca$rotation[1:100, 2:3], labels=d_names$last_name[1:100], cex=0.5)
plot(pca$sdev, type='l')

idx = which(d_names$last_name == "HERNANDEZ")
dists = rowSums(abs(ngrams %r-% ngrams[idx, ]))
d_names$last_name[order(dists, decreasing=FALSE)] |> head(100)

# cl = cluster::clara(ngrams, k=5, metric="manhattan", cluster.only=TRUE)

grp = with(d_names, split(last_name, soundex))
grp = grp[order(lengths(grp), decreasing=TRUE)]


d = pseudo_vf |>
    mutate(soundex = stringdist::phonetic(last_name, useBytes=TRUE),
           fl = str_sub(last_name, 1, 1))
p_r = with(d, prop.table(table(race)))

rpr = bisg(~ nm(last_name) + zip(zip), data=d, p_r=p_r)
p_xr = with(d, prop.table(table(turnout, race), 2)) |> round(3)

res = list(
    m0 = birdie(rpr, turnout ~ 1, data=d),
    m1_zip = birdie(rpr, turnout ~ proc_zip(zip), data=d),
    m1_sdx = birdie(rpr, turnout ~ soundex, data=d),
    m2_zip = birdie(rpr, turnout ~ (1 | proc_zip(zip)), data=d),
    m2_sdx = birdie(rpr, turnout ~ (1 | soundex), data=d)
    # m2_both = birdie(rpr, turnout ~ soundex + (1 | proc_zip(zip)), data=d)
)
res$m2_both =  birdie(rpr, turnout ~ fl + (1 | proc_zip(zip)), data=d)

unlist(lapply(res, function(m) {
    colSums(abs(coef(m) - p_xr)*0.5) %*% p_r
}))

do.call(rbind, lapply(res, function(m) {
    coef(m, subgroup=T) |> apply(1:2, sd) |> colMeans()
}))
