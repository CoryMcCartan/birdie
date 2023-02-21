library(here)
library(arrow)
library(tidyverse)
library(censable)

# chunk CSV
path = "data-raw/VM2--FL--2022-12-31/VM2--FL--2022-12-31-DEMOGRAPHIC.tab"
ds = open_dataset(path, format="tsv")

ds |>
    select(Voters_FirstName,
           Voters_LastName,
           Voters_FIPS,
           Residence_Addresses_CensusTract,
           Residence_Addresses_CensusBlock,
           Voters_Gender,
           CountyEthnic_Description,
           Parties_Description) |>
    group_by(Voters_FIPS) |>
    write_dataset("data-raw/VM2--FL--2022-12-31/parquet/", format="parquet")

d = open_dataset("data-raw/VM2--FL--2022-12-31/parquet/") |>
    filter(Voters_FIPS == 86, CountyEthnic_Description != "") |> # Miami-Dade
    collect() |>
    mutate(across(matches("Voters_.*Name"), ~ iconv(.x,"UTF-8", "UTF-8", sub = ""))) |>
    transmute(first_name = Voters_FirstName,
              last_name = Voters_LastName,
              county = str_c(match_fips("FL"), str_pad(Voters_FIPS, 3, pad="0")),
              tract = Residence_Addresses_CensusTract,
              block = Residence_Addresses_CensusBlock,
              gender = Voters_Gender,
              race = recode_factor(CountyEthnic_Description,
                                   "White Self Reported" = "white",
                                   "African or Af-Am Self Reported"="black",
                                   "Hispanic" = "hisp",
                                   "East Asian"="asian",
                                   "Korean"= "asian",
                                   "Other Undefined Race" = "other",
                                   "Native American (self reported)" = "aian"),
              party = recode_factor(Parties_Description,
                                    "Democratic" = "dem",
                                    .default = "ind",
                                    "Republican" = "rep",
                                    "Libertarian" = "lib")) |>
    mutate(across(everything(), as.factor),
           GEOID = str_c(county, tract, block)) |>
    # left_join(baf, by="GEOID") |>
    select(-GEOID)

write_rds(d, "data-raw/l2_miami.rds", compress="xz")
