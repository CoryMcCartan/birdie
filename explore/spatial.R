library(raceproxy)
library(tidyverse)
library(here)
library(wacolors)

if (file.exists(voterfile <- here("data-raw/nc_voters.rds"))) {
    voters = readRDS(voterfile)
} else {
    voters = make_nc_statewide(voterfile)
}

geom = tigris::zctas(state="NC", year=2010) %>%
    rmapshaper::ms_simplify(keep_shapes=TRUE)
cities = sf::st_centroid(filter(geom, ZCTA5CE10 %in% c("27601", "28202", "27401", "28403")))

pl_raw = pl_read(pl_url("NC", 2020))
pl = pl_raw %>%
    pl_subset(sumlev="150") %>%
    pl_select_standard()
rm(pl_raw)
geom_bg = tigris::block_groups(state="NC", year=2020)
idx_bg = geomander::geo_match(geom_bg, geom, method="center")
d_cens_20 = pl %>%
    mutate(pop_asian = pop_asian + pop_nhpi,
           vap_asian = vap_asian + vap_nhpi,
           pop_other = pop_other + pop_two,
           vap_other = vap_other + vap_two) %>%
    select(-pop_nhpi, -vap_nhpi, -pop_two, -vap_two) %>%
    select(pop:vap_other) %>%
    mutate(zip = geom$ZCTA5CE10[idx_bg]) %>%
    group_by(zip) %>%
    summarize(across(everything(), sum))
d_cens_10 = censable::build_dec("zcta", "NC", year=2010, geometry=F) %>%
    mutate(zip = str_sub(GEOID, 3),
           pop_asian = pop_asian + pop_nhpi,
           vap_asian = vap_asian + vap_nhpi,
           pop_other = pop_other + pop_two,
           vap_other = vap_other + vap_two) %>%
    select(-pop_nhpi, -vap_nhpi, -pop_two, -vap_two) %>%
    select(zip, pop:vap_other)

shared_zip = intersect(d_cens_10$zip, d_cens_20$zip)

m_10 = (select(d_cens_10, starts_with("vap_")) / d_cens_10$vap)[,c(1,2,3,5,4,6)]
m_20 = (select(d_cens_20, starts_with("vap_")) / d_cens_20$vap)[,c(2,3,1,5,4,6)]
d_diff = m_20[d_cens_20$zip %in% shared_zip, ] - m_10[d_cens_10$zip %in% shared_zip, ]
names(d_diff) = str_replace(names(d_diff), "vap_", "chg_")
d_diff = as_tibble(cbind(zip=shared_zip, d_diff))

m_10 = filter(d_cens_10, zip %in% shared_zip) %>%
    select(vap_white, vap_black, vap_hisp, vap_asian, vap_aian, vap_other) %>%
    mutate(across(everything(), ~ . / sum(.)))
m_20 = filter(d_cens_20, zip %in% shared_zip) %>%
    select(vap_white, vap_black, vap_hisp, vap_asian, vap_aian, vap_other) %>%
    mutate(across(everything(), ~ . / sum(.)))
d_diff = m_20 - m_10
names(d_diff) = str_replace(names(d_diff), "vap_", "chg_")
d_diff = as_tibble(cbind(zip=shared_zip, d_diff))

d = slice_sample(voters, n=50e3) %>%
    mutate(gender = as.factor(coalesce(gender, "U"))) %>%
    filter(!is.na(age))

p_xr = prop.table(table(d$party, d$race))
p_x = rowSums(p_xr)
p_r = colSums(p_xr)

bisg = predict_race_sgz(last_name, zip, data=d, p_r=p_r)

x = numeric(nrow(bisg))
for (i in seq_along(p_r)) {
    idx = as.integer(d$race) == i
    x[idx] = bisg[[i]][idx]
}

d = d %>%
    bind_cols(bisg) %>%
    mutate(pr_race = x)


d_cal = d %>%
    group_by(zip) %>%
    summarize(white = mean(race == "white") - mean(pr_white),
              black = mean(race == "black") - mean(pr_black),
              hisp = mean(race == "hisp") - mean(pr_hisp),
              asian = mean(race == "asian") - mean(pr_asian),
              aian = mean(race == "aian") - mean(pr_aian),
              other = mean(race == "other") - mean(pr_other)) %>%
    pivot_longer(-zip, names_to="race", values_to="bias") %>%
    left_join(count(d, zip, race), by=c("zip", "race")) %>%
    left_join(summarize(group_by(d, zip), dem=mean(party == "dem"), pop=n()), by="zip") %>%
    left_join(d_diff, by="zip") %>%
    mutate(n = coalesce(n, 0),
           bias = -bias) %>%
    left_join(select(geom, zip=ZCTA5CE10, geometry), by="zip")

ggplot(d_cal, aes(fill=bias, geometry=geometry)) +
    facet_wrap(~ race, ncol=2) +
    geom_sf(size=0) +
    geom_sf(data=cities, inherit.aes=F, shape=1, size=8) +
    labs(fill="Calibration-in-the-large:\nE[r\U0302 | G] - E[R | G]") +
    scale_fill_wa_c("lopez", midpoint=0, limits=c(-0.1, 0.1),
                    oob=scales::squish, labels=\(x) scales::percent(x, 1)) +
    theme_void(base_size=16) +
    theme(legend.position="bottom",
          legend.key.width=unit(1, "cm"))

d_cal %>%
    group_by(zip) %>%
    mutate(chg = c_across(chg_white:chg_other)[c(1, 7, 13, 19, 25, 31)]) %>%
    ungroup() %>%
    select(-starts_with("chg_")) %>%
    drop_na(chg) %>%
ggplot(aes(fill=chg, geometry=geometry)) +
    facet_wrap(~ race, ncol=2) +
    geom_sf(size=0) +
    geom_sf(data=cities, inherit.aes=F, shape=1, size=8) +
    scale_fill_wa_c("lopez", midpoint=0, #limits=c(-0.1, 0.1),
                    oob=scales::squish, labels=\(x) scales::percent(x, 1)) +
    theme_void(base_size=16) +
    theme(legend.position="bottom",
          legend.key.width=unit(1, "cm"))

d_cal %>%
    mutate(area = geos::geos_area(geometry),
           dens = pop / area) %>%
    filter(pop > 0, area > 0) %>%
    group_by(zip) %>%
    mutate(chg = c_across(chg_white:chg_other)[c(1, 7, 13, 19, 25, 31)]) %>%
    ungroup() %>%
    select(-starts_with("chg_")) %>%
    drop_na(chg) %>%
ggplot(aes(chg, bias, size=n, color=race)) +
    geom_point(alpha=0.2) +
    geom_smooth(formula=y~x, method="loess", span=0.8, se=F) +
    #coord_cartesian(ylim=c(-0.1, 0.1), xlim=c(-0.01, 0.01)) +
    coord_cartesian(ylim=c(-0.15, 0.15)) +
    scale_size_area(guide="none")
