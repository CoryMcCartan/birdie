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
                                 col_types=cols(birth_age="i", birth_year="i", .default="c"))

    race_codes = c(A="asian", B="black", I="other", M="other", O="other",
                   P="other", W="white")
    party_codes = c(UNA="ind", DEM="dem", REP="rep", LIB="lib")

    voters = voters_raw %>%
        filter(race_code != "U", ethnic_code != "UN") %>%
        mutate(race = factor(if_else(ethnic_code == "HL", "hisp",
                                     race_codes[race_code]),
                             levels=c("white", "black", "hisp", "asian", "other")),
               gender = as_factor(gender_code),
               party = factor(party_codes[party_cd],
                              levels=c("dem", "ind", "rep", "lib")),
               lic = drivers_lic == "Y") %>%
        select(last_name:middle_name, suffix=name_suffix_lbl, zip=zip_code,
               race, gender, party, lic) %>%
        drop_na(last_name)

    unlink(zipfile)
    unlink(rawfile)

    voters
}


calc_joints = function(p_xr, voters, fit, warmup = 1:100) {
    xr = list(true = p_xr)

    xr$base_surn = map(rownames(p_xr), ~ colMeans(fit$base_surn * (voters$party == .))) %>%
        do.call(rbind, .) %>%
        `rownames<-`(rownames(p_xr))

    xr$base = map(rownames(p_xr), ~ colMeans(fit$baseline * (voters$party == .))) %>%
        do.call(rbind, .) %>%
        `rownames<-`(rownames(p_xr))

    if ("gibbs" %in% names(fit)) {
        xr$gibbs = map(1:5, function(race) {
            map_dbl(rownames(p_xr), ~ mean((fit$gibbs[, -warmup] == race) * (voters$party == .)))
        }) %>%
            do.call(cbind, .) %>%
            `rownames<-`(rownames(p_xr)) %>%
            `colnames<-`(colnames(p_xr))
    }

    if ("stan" %in% names(fit)) {
        if (class(fit$stan) == "stanfit") {
            draws = rstan::extract(fit$stan, "p_xr")$p_xr
            xr$stan = (colMeans(draws) %*% diag(p_r)) %>%
                `rownames<-`(rownames(p_xr)) %>%
                `colnames<-`(colnames(p_xr))
            xr$stan_low = apply(draws, 2:3, \(x) quantile(x, 0.05)) %*% diag(p_r)
            xr$stan_high = apply(draws, 2:3, \(x) quantile(x, 0.95)) %*% diag(p_r)
            coverage = mean((p_xr > xr$stan_low) & (p_xr < xr$stan_high))
            cli_inform("Average width: {round(mean(xr$stan_high - xr$stan_low), 3)}")
            cli_inform("Empirical coverage: {scales::percent(coverage)}")
        } else {
            xr$stan =  matrix(nrow=length(p_x), ncol=length(p_r))
            for (j in 1:length(p_r)) {
                for (k in 1:length(p_x)) {
                    xr$stan[k, j] = fit$stan$par[str_glue("p_xr[{k},{j}]")]
                }
            }
            xr$stan = (xr$stan %*% diag(p_r)) %>%
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
