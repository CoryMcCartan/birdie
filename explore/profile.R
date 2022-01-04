library(tidyverse)
library(microbenchmark)
library(here)

devtools::load_all(here("."))

voterfile <- here("data/nc_voters.rds")
voters = read_rds(voterfile) %>%
    mutate(id=1:n(), .before=last_name) %>%
    select(-starts_with("pred."))
d = slice_sample(voters, n=1e5, replace=TRUE)

cat("*****\n")

invisible(model_race(party, last_name, zip, data=d, iter=1e2))
