# UTILITIES for studying NC


make_nc_df = function(county="Dare") {
    rlang::check_installed("tigris", "NC voter data")
    counties = tigris::fips_codes$county[tigris::fips_codes$state == "NC"]
    county_i = match(paste(county, "County"), counties)

    url = glue::glue("https://s3.amazonaws.com/dl.ncsbe.gov/data/ncvoter{county_i}.zip")
    tmp = withr::local_tempdir()
    zipfile = paste(tmp, "ncvoters.zip")
    download.file(url, zipfile)
    unzip(zipfile, exdir=tmp)
    rawfile = glue::glue("{tmp}/ncvoter{county_i}.txt")

    voters_raw = readr::read_tsv(rawfile, show_col_types=F,
                                 col_types=readr::cols(birth_age="i",
                                                       birth_year="i",
                                                       .default="c"))

    race_codes = c(A="asian", B="black", I="other", M="other", O="other",
                   P="other", W="white")
    party_codes = c(UNA="ind", DEM="dem", REP="rep", LIB="lib")

    voters = voters_raw %>%
        filter(race_code != "U", ethnic_code != "UN") %>%
        mutate(race = factor(dplyr::if_else(ethnic_code == "HL", "hisp",
                                 race_codes[race_code]),
                             levels=c("white", "black", "hisp", "asian", "other")),
               gender = as.factor(gender_code),
               party = factor(party_codes[party_cd],
                              levels=c("dem", "ind", "rep", "lib")),
               age = cut(birth_age, breaks=c(18, 20, 25, 30, 35, 40, 45, 50, 55,
                                       60, 62, 65, 67, 70, 75, 80, 85, 150),
                         right=FALSE),
               lic = drivers_lic == "Y") %>%
        select(last_name:middle_name, suffix=name_suffix_lbl, zip=zip_code,
               race, gender, age, birth_state, party, lic) %>%
        dplyr::slice(which(!is.na(.$last_name)))

    unlink(zipfile)
    unlink(rawfile)

    voters
}


calc_joints = function(p_xr, voters, fit, warmup = 1:100) {
    xr = list(true = p_xr)

    if ("bis" %in% names(fit)) {
        xr$bis = map(rownames(p_xr), ~ colMeans(fit$bis * (voters$party == .))) %>%
            do.call(rbind, .) %>%
            `rownames<-`(rownames(p_xr))
    }

    if ("bisg" %in% names(fit)) {
        xr$bisg = map(rownames(p_xr), ~ colMeans(fit$bisg * (voters$party == .))) %>%
            do.call(rbind, .) %>%
            `rownames<-`(rownames(p_xr))
    }

    if ("gibbs" %in% names(fit)) {
        xr$gibbs = map(rownames(p_xr), ~ colMeans(fit$gibbs * (voters$party == .))) %>%
            do.call(rbind, .) %>%
            `rownames<-`(rownames(p_xr))
        #xr$gibbs = map(1:5, function(race) {
        #    map_dbl(rownames(p_xr), ~ mean((fit$gibbs[, -warmup] == race) * (voters$party == .)))
        #}) %>%
        #    do.call(cbind, .) %>%
        #    `rownames<-`(rownames(p_xr)) %>%
        #    `colnames<-`(colnames(p_xr))
    }

    if ("nonparam" %in% names(fit)) {
        if (class(fit$nonparam) == "stanfit") {
            draws = rstan::extract(fit$nonparam, "p_xr")$p_xr
            xr$nonparam = (colMeans(draws) %*% diag(p_r)) %>%
                `rownames<-`(rownames(p_xr)) %>%
                `colnames<-`(colnames(p_xr))
            xr$nonparam_low = apply(draws, 2:3, \(x) quantile(x, 0.05)) %*% diag(p_r)
            xr$nonparam_high = apply(draws, 2:3, \(x) quantile(x, 0.95)) %*% diag(p_r)
            coverage = mean((p_xr > xr$nonparam_low) & (p_xr < xr$nonparam_high))
            cli_inform("Average width: {round(mean(xr$nonparam_high - xr$nonparam_low), 3)}")
            cli_inform("Empirical coverage: {scales::percent(coverage)}")
        } else {
            xr$nonparam =  matrix(nrow=length(p_x), ncol=length(p_r))
            for (j in 1:length(p_r)) {
                for (k in 1:length(p_x)) {
                    xr$nonparam[k, j] = fit$nonparam$par[str_glue("p_xr[{k},{j}]")]
                }
            }
            xr$nonparam = (xr$nonparam %*% diag(p_r)) %>%
                `rownames<-`(rownames(p_xr)) %>%
                `colnames<-`(colnames(p_xr))
        }
    }

    if ("additive" %in% names(fit)) {
        if (class(fit$additive) == "stanfit") {
            draws = rstan::extract(fit$additive, "p_xr")$p_xr
            xr$additive = (colMeans(draws) %*% diag(p_r)) %>%
                `rownames<-`(rownames(p_xr)) %>%
                `colnames<-`(colnames(p_xr))
            xr$additive_low = apply(draws, 2:3, \(x) quantile(x, 0.05)) %*% diag(p_r)
            xr$additive_high = apply(draws, 2:3, \(x) quantile(x, 0.95)) %*% diag(p_r)
            coverage = mean((p_xr > xr$additive_low) & (p_xr < xr$additive_high))
            cli_inform("Average width: {round(mean(xr$additive_high - xr$additive_low), 3)}")
            cli_inform("Empirical coverage: {scales::percent(coverage)}")
        } else {
            xr$additive =  matrix(nrow=length(p_x), ncol=length(p_r))
            for (j in 1:length(p_r)) {
                for (k in 1:length(p_x)) {
                    xr$additive[k, j] = fit$additive$par[str_glue("p_xr[{k},{j}]")]
                }
            }
            xr$additive = (xr$additive %*% diag(p_r)) %>%
                `rownames<-`(rownames(p_xr)) %>%
                `colnames<-`(colnames(p_xr))
        }
    }

    xr
}

print_cond = function(x, title=NULL) {
    p_r = parent.frame()$p_r
    out = 100 * x %*% diag(1 / p_r)
    colnames(out) = names(p_r)
    if (!is.null(title)) cat(toupper(title), "\n")
    print(round(out))
}

plot_county_puma = function(county=NULL) {
    suppressMessages({
    library(ggplot2)
    d_puma = tigris::pumas("NC", cb=TRUE)
    d_cty = tigris::counties("NC", cb=TRUE)
    d_zip = tigris::zctas(cb=F, state="NC", year=2010)
    })
    names = abbreviate(d_cty$NAME, 5, dot=T)

    if (!is.null(county))
        bbox = sf::st_bbox(filter(d_cty, NAME %in% county))
    else
        bbox = sf::st_bbox(d_cty)

    p = ggplot(d_cty) +
        geom_sf(size=0.75, color="black") +
        geom_sf(data=d_zip, size=0.1, fill=NA, color="yellow") +
        geom_sf(data=d_puma, size=0.25, fill=NA, color="red") +
        geom_sf_text(aes(label=names), size=2.5) +
        coord_sf(xlim=bbox[c(1,3)], ylim=bbox[c(2,4)]) +
        theme_void()
    suppressWarnings(print(p))
}


prep_age_sex_race_zip = function() {
    sf1v = tidycensus::load_variables(2010, dataset="sf1") %>%
        rename(variable=name)
    tables = c("P012", "P012I", "P012B", "P012D", "P012H")
    d_raw = map_dfr(tables, ~ tidycensus::get_decennial("zcta", table=., state="NC"))
    d = d_raw %>%
        select(-NAME) %>%
        left_join(sf1v, by="variable") %>%
        filter(!label %in% c("Total", "Total!!Female", "Total!!Male")) %>%
        mutate(zip = str_sub(GEOID, 3),
               race = case_when(str_detect(concept, "WHITE") ~ "white",
                                str_detect(concept, "BLACK") ~ "black",
                                str_detect(concept, "ASIAN") ~ "asian",
                                str_detect(concept, "\\(HISPANIC") ~ "hisp",
                                TRUE ~ "total")) %>%
        select(-concept, -GEOID) %>%
        pivot_wider(id_cols=c(zip, label), names_from=race, values_from=value) %>%
        separate(label, c(NA, "sex", "age1", "age2"), sep="(!!|( to | and ))", fill="right") %>%
        mutate(age1 = parse_number(age1),
               age2 = coalesce(parse_number(age2), 149),
               age2 = if_else(age1 %in% 20:21, 24, age2),
               age1 = if_else(age1 %in% 21:22, 20, age1),
               age = str_glue("[{age1},{age2+1L})")) %>%
        suppressWarnings() %>%
        filter(age1 >= 18) %>%
        mutate(other = total - white - black - asian - hisp,
               sex = str_sub(sex, 1, 1)) %>%
        select(zip, age, sex, total, white, black, hisp, asian, other) %>%
        group_by(zip, age, sex) %>%
        summarize(across(total:other, sum)) %>%
        ungroup()
    write_rds(d, "data/nc_zip_age_sex_race.rds", compress="xz")
}


# Helper function to make an R|G table
census_zip_age_sex_table = function(GZ, GZ_vec, p_r, regularize=TRUE, count=FALSE) {
    if (!rlang::is_installed("zipWRUext2"))
        cli_abort(c("{.pkg zipWRUext2} must be installed to use ZIP code information automatically.",
                    ">"=' {.code devtools::install_github("https://github.com/jcuriel-unc/zipWRUext",
                    subdir="zipWRUext2")} to install.'))
    if (length(p_r) != 5)
        cli_abort(c("Number of racial categories doesn't match the Census data.",
                    "i"="Categories should be White, Black, Hispanic, Asian, and other."))

    x = read_rds("data/nc_zip_age_sex_race.rds")
    x_nozip = group_by(x, age, sex) %>%
        summarize(across(total:other, sum)) %>%
        ungroup() %>%
        mutate(zip="<none>")
    x = bind_rows(x, x_nozip)
    x_nosex = group_by(x, zip, age) %>%
        summarize(across(total:other, sum)) %>%
        ungroup() %>%
        mutate(sex="U")
    x = bind_rows(x, x_nosex)

    GZ$order = as.integer(GZ_vec)
    GZ$gender = as.character(GZ$gender)
    GZ$zip[!GZ$zip %in% unique(x$zip)] = "<none>"
    GZ = distinct(GZ, zip, age, gender, order)
    x = inner_join(x, GZ, by=c("zip", "age", "sex"="gender")) %>%
        arrange(order) %>%
        select(-order)

    alpha = if (regularize) c(3.1, 0.6, 0.8, 0.3, 0.2) else rep(0.001, 5)
    for (i in seq_along(p_r)) {
        denom = if (count) 1 else (x[[4]] + sum(alpha))
        x[[4+i]] = (if_else(x[[4+i]] < 0, 0, x[[4+i]]) + alpha[i]) / denom
    }

    x[, -4]
}
