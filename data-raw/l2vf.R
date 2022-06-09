library(tidyverse)
library(censable)
library(here)

states <- c("AL", "FL", "GA", "LA", "NC", "SC")

# prepare ZCTA BAF ------
url_baf <- "https://www2.census.gov/geo/docs/maps-data/data/rel2020/zcta520/tab20_zcta520_tabblock20_natl.txt"
path_baf <- "data-raw/zcta_baf.txt"
if (!file.exists(here(path_baf))) download.file(url_baf, here(path_baf))

regex_fips <- str_c("(", paste(map_chr(states, match_fips), collapse="|"), ")")
baf <- read_delim(path_baf, delim="|", col_types="_c_______c_______") %>%
    rename(zip=GEOID_ZCTA5_20, GEOID=GEOID_TABBLOCK_20) %>%
    drop_na(zip, GEOID) %>%
    mutate(zip=as.factor(zip)) %>%
    filter(str_starts(GEOID, regex_fips))

# process voter files ------
# input: unzipped L2 voter files in data-raw/VF
for (state in states) {
    cat("Processing", state, "\n")
    path_raw = here(str_glue("data-raw/VF/{state}.csv"))
    path_tidy = here(str_glue("data-raw/VF/l2vf_{tolower(state)}.rds"))
    read_csv(path_raw, col_types=cols(.default="c")) %>%
        mutate(across(matches("Voters_.*Name"), ~ iconv(.x,"UTF-8", "UTF-8", sub = ""))) %>%
        drop_na(CountyEthnic_Description) %>%
        transmute(first_name = Voters_FirstName,
                  last_name = Voters_LastName,
                  county = str_c(match_fips(state), Voters_FIPS),
                  tract = Residence_Addresses_CensusTract,
                  block = Residence_Addresses_CensusBlock,
                  gender = Voters_Gender,
                  race = recode(CountyEthnic_Description,
                                "White Self Reported" = "white",
                                "African or Af-Am Self Reported"="black",
                                "Hispanic" = "hisp",
                                "East Asian"="asian",
                                "Korean"= "asian",
                                "Other Undefined Race" = "other",
                                "Native American (self reported)" = "aian")) %>%
        mutate(across(everything(), as.factor),
               GEOID = str_c(county, tract, block)) %>%
        left_join(baf, by="GEOID") %>%
        select(-GEOID) %>%
        write_rds(path_tidy, compress="gz")
    file.remove(path_raw)
}

file.remove(path_baf)
