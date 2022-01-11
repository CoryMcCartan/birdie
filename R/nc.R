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
        mutate(race = as_factor(if_else(ethnic_code == "HL", "hisp",
                                        race_codes[race_code])),
               gender = as_factor(gender_code),
               party = as_factor(party_codes[party_cd]),
               lic = drivers_lic == "Y") %>%
        select(last_name:middle_name, suffix=name_suffix_lbl, zip=zip_code,
               race, gender, party, lic) %>%
        drop_na(last_name)

    unlink(zipfile)
    unlink(rawfile)

    voters
}
