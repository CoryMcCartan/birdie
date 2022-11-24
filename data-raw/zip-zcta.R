# Prepare ZIP-ZCTA crosswalk
library(easycensus)
library(tidyverse)
library(here)

zip_raw = read_csv("data-raw/zip_zcta_xref.csv", show_col_types=FALSE)

zip_xw = zip_raw |>
    select(zip=zip_code, zcta) |>
    arrange(zip, zcta)

# run sysdata.R to save!

