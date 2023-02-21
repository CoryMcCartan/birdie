library(dplyr)
library(forcats)
library(readr)
library(stringr)
library(stringdist)
library(here)
devtools::load_all(".")

if (!file.exists(path <- here("data-raw/ipums-1930/proc_light.rds"))) {
    library(ipumsr)
    ddi <- read_ipums_ddi(here("data-raw/ipums-1930/usa_00002.xml"))
    raw <- read_ipums_micro(ddi)

    d <- raw |>
        select(wt=PERWT, name=NAMELAST, race=RACE, race_d=RACED,
               hisp=HISPAN, hisp_d=HISPAND, hisp_sur=SPANNAME, tribe=TRIBE,
               bpl=BPL, bpl_m=MBPL, bpl_f=FBPL, yr_img=YRIMMIG) |>
        mutate(across(race:bpl_f, haven::as_factor),
               yr_img = na_if(as.integer(yr_img), 0L)) |>
        filter(name != "", !str_detect(name, "[-=\\[\\]]"))

    nms = unique(d$name)
    nms_proc = proc_name(nms)
    d$name = nms_proc[match(d$name, nms)]
    d <- d |>
        filter(name != "", !is.na(name)) |>
        group_by(across(-wt)) |>
        summarize(wt = sum(wt), .groups="drop")

    write_rds(d, path, compress="gz")
} else {
    d = read_rds(path)
}

extra_asian <- read_csv(here("data-raw/asian_surnames.csv"), show_col_types=FALSE) |>
    mutate(name = str_to_upper(name))

cens <- readRDS(system.file("extdata", "names_2010_counts.rds",
                            package="birdie", mustWork=TRUE)) |>
    pull(last_name)
if (!file.exists(path <- here("data-raw/ipums-1930/cens_matched.rds"))) {
    idx_matched = cens %in% c(d$name, extra_asian$name)
    matched = cens[idx_matched]

    repl = matched[amatch(cens[!idx_matched], matched,
                          useBytes=TRUE, maxDist=Inf, method="lcs")]
    dists = stringdist(cens[!idx_matched], repl, method="lv")
    repl[dists > 6] = NA_character_
    cens_matched = cens
    cens_matched[!idx_matched] = repl
    write_rds(cens_matched, path, compress="gz")
} else {
    cens_matched = read_rds(path)
}

dom_bpl = c(state.name, "United States, ns", "District of Columbia",
            "Native American", "Abroad (unknown) or at sea")

d <- d |>
    mutate(
        race = fct_recode(fct_relabel(race, str_to_lower),
            "black" = "black/african american",
            "aian" = "american indian or alaska native",
            "asian" = "other asian or pacific islander"
        ),
        hisp = fct_recode(fct_relabel(hisp, str_to_lower),
            "no"="not hispanic",
            "prq"="puerto rican",
            "hisp"="other"
        ),
        race = if_else(hisp == "no", race,
                       if_else(hisp_d == "Spanish", "spn", hisp)),
        ancestry_d = factor(case_when(
            race == "aian" ~ "Native",
            race_d == "American Indian/Alaska Native" ~ "Native",
            str_starts(race_d, "Hawaiian") ~ "asian", # other asian, pac isl.
            race == "white" & (bpl == "AFRICA" | bpl_f == "AFRICA") ~ "South Africa",
            !bpl %in% dom_bpl ~ bpl,
            !bpl_f %in% dom_bpl ~ bpl_f,
            race == "chinese" ~ "China",
            race == "japanese" ~ "Japan",
            race == "spn" ~ "Spain",
            race == "prq" ~ "Puerto Rico",
            race == "cuban" ~ "Cuba",
            race == "mexican" ~ "Mexico",
            race_d == "Mexican (1930)" ~ "Mexico",
            race == "black" ~ "Black",
            race_d == "Black/African American" ~ "Black",
            race_d == "Asian Indian (Hindu 1920_1940)" ~ "India",
            race_d == "Filipino" ~ "Philippines",
            race_d == "Portuguese" ~ "Portugal",
            race_d == "Korean" ~ "Korea",
            race_d == "Other Asian or Pacific Islander (1920,1980)" ~ "asian",
            TRUE ~ "USA"
        )),
        ancestry = fct_infreq(fct_collapse(ancestry_d,
            anglo = c("USA", "England", "Wales", "Canada", "United Kingdom, ns",
                      "Scotland", "Australia and New Zealand", "South Africa",
                      "St. Pierre and Miquelon", "Gibraltar"),
            black = "Black",
            china = "China",
            japkor = c("Japan", "Korea"),
            seasia = c("Singapore", "Malaysia", "Thailand",
                      "Southeast Asia, ns", "Indonesia"),
            sasia = c("India", "Southwest Asia, nec/ns"),
            asiapac = c("asian", "Asia, nec/ns", "Pacific Islands", "Philippines",
                        "Guam", "American Samoa"),
            native = "Native",
            mena = c("Afghanistan", "Asia Minor, n.s.", "Asia Minor, ns",
                     "Iran", "Iraq", "Syria", "Turkey", "Lebanon",
                     "Middle East, ns", "Persian Gulf States, n.s.",
                     "Yemen Arab Republic (North)"),
            afrocar = c("AFRICA", "U.S. Virgin Islands", "Atlantic Islands",
                        "US Virgin Islands", "West Indies"),
            italy = c("Italy", "San Marino", "Malta"),
            irel = "Ireland",
            weur = c("Spain", "Portugal", "France",
                     "Western Europe, ns"),
            german = c("Germany", "Netherlands", "Belgium", "Luxembourg",
                       "Switzerland", "Austria", "Liechtenstein",
                       "Europe, ns", "Europe, nec/ns", "Northern Europe, ns"),
            eeur = c("Czechoslovakia", "Poland", "Hungary", "Romania",
                     "Bulgaria", "Yugoslavia", "Central Europe, ns",
                     "Estonia", "Latvia", "Lithuania", "Eastern Europe, ns",
                     "Baltic States, ns"),
            agean = c("Greece", "Cyprus", "Albania"),
            mexico = "Mexico",
            cuba = "Cuba",
            prq = "Puerto Rico",
            latin = c("SOUTH AMERICA", "Central America"),
            nordic = c("Denmark", "Iceland", "Norway", "Sweden", "Finland", "Lapland, n.s."),
            jwsrus = c("Other USSR/Russia", "Israel/Palestine")
        ))
    ) |>
    select(name, wt, race, hisp, bpl:ancestry)

prior_ancs = prop.table(table(d$ancestry))

d_grp = collapse::rsplit(d, d$name) |>
    lapply(function(x) {
        tab = collapse::qtab(x$ancestry, w=x$wt)
        rel = prop.table(tab) / prior_ancs
        names(which.max(rel))
    }) |>
    do.call(rbind, args=_) |>
    as.data.frame() |>
    tibble::rownames_to_column("name") |>
    as_tibble() |>
    rename(group=V1)

d_grp_extra = extra_asian |>
    left_join(d_grp, by="name", suffix=c("", "_30")) |>
    mutate(group = fct_collapse(
        group,
        japkor=c("japan", "korea"),
        sasia="india",
        asiapac="philippines",
        seasia="vietnam"
        )) |>
    filter(!(group == "asiapac" & group_30 %in% c("mexico", "cuba", "prq"))) |>
    select(-group_30)

d_grp = d_grp |>
    anti_join(d_grp_extra, by="name") |>
    bind_rows(d_grp_extra) |>
    mutate(group = factor(group, levels=names(prior_ancs))) |>
    arrange(name)

d_grp = tibble(last_name=cens, name=cens_matched) |>
    full_join(d_grp, by="name") |>
    mutate(last_name = coalesce(last_name, name)) |>
    select(-name)


d_surname = collapse::rsplit(d, d$name) |>
    lapply(function(x) {
        rel = prop.table(collapse::qtab(x$ancestry, w=x$wt)) / prior_ancs
        names(which.max(rel))
    }) |>
    do.call(rbind, args=_) |>
    as.data.frame() |>
    tibble::rownames_to_column("name") |>
    as_tibble()

d_surname = collapse::rsplit(d, d$name) |>
    lapply(function(x) {
        prop.table(collapse::qtab(x$ancestry, w=x$wt))
    }) |>
    do.call(rbind, args=_) |>
    as.data.frame() |>
    tibble::rownames_to_column("name") |>
    as_tibble()

d_surname = tibble(last_name=cens, name=cens_matched) |>
    full_join(d_surname, by="name") |>
    mutate(last_name = coalesce(last_name, name)) |>
    select(-name) |>
    mutate(across(-last_name, ~ coalesce(., 0) |> round(2)))
write_rds(d_surname, here("data-out/surn_cov.rds"), compress="xz")


# exploratory
top_surn <- function(col, val, n=50) {
    x = rlang::eval_tidy(rlang::enquo(col), d)
    nms = d$name[as.integer(x) == pmatch(val, levels(x), duplicates.ok=TRUE)]
    table(nms) |>
        sort(decreasing=TRUE) |>
        head(n) |>
        names()
}

sort(table(as.character(d$bpl_f)), decreasing=T)
top_surn(bpl_f, "Afghanistan")

nms = d$name[d$race=="Japanese" & (d$bpl == "Japan" | d$bpl_f == "Japan")] |> unique()


voters = voters |>
    left_join(d_surname, by="last_name") |>
    mutate(across(anglo:sasia, ~ coalesce(., 0)))

d_cor_dem = voters |>
    filter(!is.na(party)) |>
    group_by(race) |>
    summarize(across(anglo:sasia, ~ cov(., party =="dem"))) |>
    tidyr::pivot_longer(-race, names_to="group")
d_cor_rep = voters |>
    filter(!is.na(party)) |>
    group_by(race) |>
    summarize(across(anglo:sasia, ~ cov(., party =="rep"))) |>
    tidyr::pivot_longer(-race, names_to="group")
d_cor_turn = voters |>
    group_by(race) |>
    summarize(across(anglo:sasia, ~ cov(., turnout))) |>
    tidyr::pivot_longer(-race, names_to="group")
d_cor_lic = voters |>
    group_by(race) |>
    summarize(across(anglo:sasia, ~ cov(., lic))) |>
    tidyr::pivot_longer(-race, names_to="group")

library(ggplot2)
library(wacolors)
ggplot(d_cor_dem, aes(race, fct_rev(fct_inorder(group)), fill=value)) +
    geom_tile() +
    scale_fill_wa_c("vantage", midpoint=0) +
    labs(x="Race", y="Surname-based group", fill="Covariance",
         title="Within-race covariance with Dem. registration")

