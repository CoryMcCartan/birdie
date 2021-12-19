library(tidyverse)
library(here)
library(zipWRUext2)

if (!file.exists(voterfile <- here("data/nc_voters.rds"))) {
    url = "https://s3.amazonaws.com/dl.ncsbe.gov/data/ncvoter31.zip"
    zipfile = "data/ncvoter31.zip"
    download.file(url, here(zipfile))
    unzip(here(zipfile), exdir=here("data"))
    rawfile = here("data/ncvoter31.txt")

    voters_raw = read_tsv(rawfile, show_col_types=F,
                          col_types=cols(birth_age="i", birth_year="i", .default="c"))

    race_codes = c(A="asian", B="black", I="other", M="other", O="other",
                   P="other", W="white")
    party_codes = c(UNA="ind", DEM="dem", REP="rep", LIB="lib")

    voters = voters_raw %>%
        filter(race_code != "U", ethnic_code != "UN") %>%
        mutate(race = as_factor(if_else(ethnic_code == "HL", "hisp",
                                        race_codes[race_code])),
               gender = as_factor(gender_code),
               party = as_factor(party_codes[party_cd]),
               lic = drivers_lic == "Y") %>%
        select(last_name:middle_name, suffix=name_suffix_lbl, zip=zip_code,
               race, gender, party, lic)

    write_rds(voters, voterfile, compress="gz")
    unlink(zipfile)
    unlink(rawfile)
} else {
    voters = read_rds(voterfile)
}

voters = mutate(voters, id=1:n(), .before=last_name)

voters_imp = zip_wru(voters, "NORTH CAROLINA", year1=2018,
                     zip_col="zip", surname_field="last_name") %>%
    as_tibble()

voters = left_join(voters, select(voters_imp, id, starts_with("pred.")), by="id")

P_r = c(p_whi=0.626, p_bla=0.222, p_his=0.098, p_asi=0.032, p_oth=0.022)
surnames = as.character(wru::surnames2010$surname)
surn_used = surnames %in% voters$last_name
P_rs = t(as.matrix(wru::surnames2010[surn_used, -1])) %>%
    cbind(P_r)
surnames = c(surnames[surn_used], "<NA>")
rownames(P_rs) = c("white", "black", "hisp", "asian", "other")

voters = voters %>%
    mutate(surname=factor(if_else(last_name %in% surnames, last_name, "<NA>"),
                          level=surnames))
P_s = prop.table(with(voters, table(surname)))
P_rs = P_rs %*% diag(P_s)
colnames(P_rs) = surnames
P_xs = prop.table(with(voters, table(party, surname)))

P_xr = prop.table(with(voters, table(party, race)))
P_xr = P_xr[, match(rownames(P_rs), colnames(P_xr))]
P_xr_low = outer(rowSums(P_xr), colSums(P_xr), \(x, y) pmax(0, x+y-1))
P_xr_high = outer(rowSums(P_xr), colSums(P_xr), pmin)

P_xr_est1 = P_xs %*% MASS::ginv(P_rs) %*% diag(rowSums(P_rs))
P_xr_est1[P_xr_est1 > P_xr_high] = P_xr_high[P_xr_est1 > P_xr_high]
P_xr_est1[P_xr_est1 < P_xr_low] = P_xr_low[P_xr_est1 < P_xr_low]
# rake
for (i in 1:5) {
    P_xr_est1 = P_xr_est1 %*% diag(colSums(P_xr) / colSums(P_xr_est1))
    P_xr_est1 = diag(rowSums(P_xr) / rowSums(P_xr_est1)) %*% P_xr_est1
}

P_xr_est2 = matrix(nrow=nrow(P_xr), ncol=ncol(P_xr), dimnames=dimnames(P_xr))
for (party in rownames(P_xr)) {
    for (race in colnames(P_xr)) {
        P_xr_est2[party, race] = weighted.mean(voters$party == party, P_rs[race, voters$surname])
    }
}
P_xr_est2 = P_xr_est2 %*% diag(rowSums(P_rs))
P_xr_est2[P_xr_est2 > P_xr_high] = P_xr_high[P_xr_est2 > P_xr_high]
P_xr_est2[P_xr_est2 < P_xr_low] = P_xr_low[P_xr_est2 < P_xr_low]
# rake
for (i in 1:5) {
    P_xr_est2 = P_xr_est2 %*% diag(colSums(P_xr) / colSums(P_xr_est2))
    P_xr_est2 = diag(rowSums(P_xr) / rowSums(P_xr_est2)) %*% P_xr_est2
}


sqrt(mean((P_xr_est1 - P_xr)^2))
sqrt(mean((P_xr_est2 - P_xr)^2))

# calibration
if (F) {
voters %>%
    mutate(prob = round(pred.bla * 20)/20) %>%
    group_by(prob) %>%
    summarize(act = mean(race == "black"),
              n = n()) %>%
ggplot(aes(prob, act, size=n)) +
    geom_abline(slope=1, color="red") +
    geom_point()
}

p(X, S) = sum[R] p(X, S, R) = sum[R] p(X | R, S) p(S | R) p(R) = sum[R] p(X | R) p(R, S)
