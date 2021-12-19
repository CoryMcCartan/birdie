library(tidyverse)
library(here)
library(wacolors)

if (file.exists(voterfile <- here("data-raw/nc_voters.rds"))) {
    voters = readRDS(voterfile)
} else {
    voters = make_nc_statewide(voterfile)
}

geom = tigris::zctas(year=2010, state="NC") %>%
    rmapshaper::ms_simplify(keep_shapes=T)
d = voters %>%
    mutate(is_irish = str_starts(last_name, "(MC|MAC)") |
               last_name %in% c("OLEARY", "OBRIEN", "CALLAGHAN", "OCONNOR", "ODONNELL", "OKEEFE", "ONEIL")) %>%
    group_by(zip, race) %>%
    summarize(irish = mean(is_irish)) %>%
    left_join(select(geom, zip=ZCTA5CE10, geometry), by="zip")


filter(d, race=="white") %>%
ggplot(aes(fill=irish, geometry=geometry)) +
    geom_sf(size=0) +
    scale_fill_wa_c(trans="sqrt") +
    theme_void()
