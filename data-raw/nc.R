# UTILITIES for studying NC

make_nc_statewide = function(voterfile) {
    counties = tigris::fips_codes$county[tigris::fips_codes$state == "NC"] %>%
        stringr::str_sub(end=-8)

    voters = bind_rows(lapply(counties, function(county) {
        cat(county, "\n")
        cty_voters = make_nc_df(county)
        cty_vhist = make_nc_vhist(county)

        left_join(cty_voters, cty_vhist, by="regnum") %>%
            mutate(across(starts_with("voted_"), ~ coalesce(., FALSE)),
                   n_voted = as.integer(rowSums(across(starts_with("voted_")))),
                   .after="lic") %>%
            select(-middle_name, -suffix, -address, -birth_state)
    })) %>%
        mutate(across(starts_with("voted_"), ~ coalesce(., FALSE)))

    voters = voters %>%
        mutate(across(where(is.character), factor))
    if (!is.null(voterfile)) {
        saveRDS(voters, voterfile, compress="xz")
    }
    voters
}

# Download and format voter file from NC
make_nc_df = function(county="Dare") {
    rlang::check_installed("tigris", "NC voter data")
    counties = tigris::fips_codes$county[tigris::fips_codes$state == "NC"]
    county_i = match(paste(county, "County"), counties)

    url = glue::glue("https://s3.amazonaws.com/dl.ncsbe.gov/data/ncvoter{county_i}.zip")
    tmp = withr::local_tempdir()
    zipfile = paste(tmp, "ncvoters.zip")
    download.file(url, zipfile, quiet=TRUE)
    unzip(zipfile, exdir=tmp)
    rawfile = glue::glue("{tmp}/ncvoter{county_i}.txt")

    voters_raw = readr::read_tsv(rawfile, show_col_types=F,
                                 col_types=readr::cols(age_at_year_end="i",
                                                       birth_year="i",
                                                       .default="c"))

    race_codes = c(A="asian", B="black", I="aian", M="other", O="other",
                   P="other", W="white")
    party_codes = c(UNA="ind", DEM="dem", REP="rep", LIB="lib")

    voters = voters_raw %>%
        filter(race_code != "U", ethnic_code != "UN") %>%
        mutate(race = factor(dplyr::if_else(ethnic_code == "HL", "hisp",
                                 race_codes[race_code]),
                             levels=c("white", "black", "hisp", "asian", "aian", "other")),
               gender = as.factor(gender_code),
               party = factor(party_codes[party_cd],
                              levels=c("dem", "ind", "rep", "lib")),
               age = cut(age_at_year_end,
                         breaks=c(18, 20, 25, 30, 35, 40, 45, 50, 55,
                                  60, 62, 65, 67, 70, 75, 80, 85, 150),
                         right=FALSE),
               address = dplyr::na_if(res_street_address, "REMOVED"),
               city = res_city_desc,
               county_name = county,
               county = paste0("37", tigris::fips_codes$county_code[county_i]),
               lic = drivers_lic == "Y") %>%
        select(regnum=voter_reg_num,
               last_name:middle_name, suffix=name_suffix_lbl,
               address, city, zip=zip_code, county, county_name,
               race, gender, age, birth_state, party, lic) %>%
        dplyr::slice(which(!is.na(.$last_name)))

    unlink(zipfile)
    unlink(rawfile)

    voters
}

make_nc_vhist = function(county = "Dare", years=2012:2022, general_only=TRUE) {
    rlang::check_installed("tigris", "NC voter data")
    counties = tigris::fips_codes$county[tigris::fips_codes$state == "NC"]
    county_i = match(paste(county, "County"), counties)

    url = glue::glue("https://s3.amazonaws.com/dl.ncsbe.gov/data/ncvhis{county_i}.zip")
    tmp = withr::local_tempdir()
    zipfile = paste(tmp, "ncvoters.zip")
    download.file(url, zipfile, quiet=TRUE)
    unzip(zipfile, exdir=tmp)
    rawfile = glue::glue("{tmp}/ncvhis{county_i}.txt")

    vhist_raw = readr::read_tsv(rawfile, show_col_types=F,
                                col_types=readr::cols(voter_reg_num="c",
                                                      election_desc="c",
                                                      .default="_"))

    months = if (general_only) "11" else formatC(1:12, width=2, flag="0")

    vhist = vhist_raw %>%
        rename(regnum=voter_reg_num) %>%
        tidyr::separate(election_desc, c("date", NA), sep=10) %>%
        tidyr::separate(date, c("month", NA, "year"), sep=c(2, 6)) %>%
        mutate(election_date = stringr::str_c("voted_", year, "_", month),
               voted = TRUE) %>%
        filter(year %in% as.character(years),
               month %in% months) %>%
        select(-month, -year) %>%
        distinct() %>%
        tidyr::pivot_wider(names_from=election_date, names_sort=TRUE,
                           values_from=voted, values_fill=FALSE)

    unlink(zipfile)
    unlink(rawfile)

    vhist
}


# AGE / SEX / RACE / ZIP TABLE

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
    if (length(p_r) != 5)
        cli_abort(c("Number of racial categories doesn't match the Census data.",
                    "i"="Categories should be White, Black, Hispanic, Asian, and other."))

    x = readr::read_rds("data-raw/nc_zip_age_sex_race.rds")
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
    x = inner_join(GZ, x, by=c("zip", "age", "gender"="sex")) %>%
        #arrange(order) %>% # probably not needed but just in case
        select(-order)

    alpha = if (regularize) c(3.1, 0.6, 0.8, 0.3, 0.2) else rep(0.001, 5)
    for (i in seq_along(p_r)) {
        denom = if (count) 1 else (x[[4]] + sum(alpha))
        x[[4+i]] = (if_else(x[[4+i]] < 0, 0, x[[4+i]]) + alpha[i]) / denom
    }

    as_tibble(x[, -4])
}
