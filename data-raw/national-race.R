#############
# Prepare national Census race time series
library(easycensus)
library(tidyverse)
library(here)
devtools::load_all(".")

years_pop = 2005:2024
get_natl_race_pop <- function(year) {
    if (year == 2010) {
        x = census_race_geo_table("us", year=2010, survey="dec", counts=FALSE)
    } else if (year == 2020) {
        return(NULL)
    } else {
        x = census_race_geo_table("us", year=year, survey="acs1", counts=FALSE)
    }

    out = unlist(x[c("white", "black", "hisp", "asian", "aian", "other")],
                 use.names=FALSE)
    names(out) = c("white", "black", "hisp", "asian", "aian", "other")
    out
}

l_race_pop = map(years_pop, get_natl_race_pop)
names(l_race_pop) = years_pop

# The standard 2020 ACS 1-year estimates were not released. Interpolate rather
# than splice a decennial estimate into this ACS series. Changes to the ACS race
# question and response coding also create a real series break between the 2019
# and 2021 estimates, most visibly in the residual/multiracial category.
l_race_pop[["2020"]] = l_race_pop[["2019"]]*0.5 + l_race_pop[["2021"]]*0.5

# Selected Population Profiles give each group's voting-age population. The
# API product begins in 2009. These population-group codes are mutually
# exclusive and match birdie's race categories.
years_vap = 2009:2024
race_popgroups = list(
    white="451",
    black="453",
    hisp="400",
    asian=c("457", "459"),
    aian="455",
    other=c("461", "463")
)

get_natl_race_vap = function(year) {
    if (year %in% c(2010, 2020))
        return(NULL)
    message("Downloading ", year, " voting-age race data")

    vars = censusapi::listCensusMetadata(
        name="acs/acs1/spp", vintage=year, type="variables"
    )
    is_vap =
        str_detect(vars$name, "^S0201_0*26E$") &
        str_detect(vars$label,
                   "^Estimate!!SEX AND AGE!!(Total population!!)?18 years and over$")
    vap_var = vars$name[is_vap]
    stopifnot(length(vap_var) == 1L)

    d = censusapi::getCensus(
        name="acs/acs1/spp", vintage=year,
        key=Sys.getenv("CENSUS_API_KEY"),
        vars=vap_var,
        region="us:1", POPGROUP="*"
    )
    get_group_vap = function(popgroups) {
        sum(d[[vap_var]][d$POPGROUP %in% popgroups])
    }
    out = map_dbl(race_popgroups, get_group_vap)
    out / sum(out)
}

l_race_vap = map(years_vap, get_natl_race_vap)
names(l_race_vap) = years_vap

# Maintain a consistent ACS series in years without standard one-year data.
l_race_vap[["2010"]] = l_race_vap[["2009"]]*0.5 + l_race_vap[["2011"]]*0.5
l_race_vap[["2020"]] = l_race_vap[["2019"]]*0.5 + l_race_vap[["2021"]]*0.5

l_race_year = map(years_pop, function(year) {
    year = as.character(year)
    list(pop=l_race_pop[[year]], vap=l_race_vap[[year]])
})
names(l_race_year) = years_pop

# Preserve the existing ZIP-ZCTA crosswalk when updating internal data.
zip_xw = get("zip_xw", asNamespace("birdie"))
message("Saving national race series to R/sysdata.rda")
save(zip_xw, l_race_year, file=here("R/sysdata.rda"), compress="xz", version=2)
