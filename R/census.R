
#' Download Census Race Data
#'
#' Downloads and prepares race-by-geography tables from U.S. census data, using
#' the [`easycensus`][easycensus::easycensus] package. Requires that an api key
#' be set up through [easycensus::cens_auth()] in that package. Supports data
#' from the decennial census and the American Community Survey at a variety of
#' levels of geographic detail. The output of this function can be used directly
#' in [bisg()].
#'
#' @param geo The geographic level to return. Common options are listed in the
#'   function signature, but any of the geographies listed at
#'   [easycensus::cens_geo()] may be used.
#' @param ... Further subgeographies to return, as in [easycensus::cens_geo()].
#' @param year The year for the data
#' @param survey The data product to use: either the decennial census (`"dec"`),
#'   or the the 1-year or 5-year ACS.
#' @param GEOIDs If `TRUE`, return the `GEOID` column as the unique geographic
#'   identifier; if `FALSE`, return a human-readable name. For example, with
#'   `geo="state"`, setting `GEOIDs=FALSE` will return a column named `state`
#'   with entries like `"Massachusetts"`.
#' @param counts If `TRUE`, return the table as actual population counts; if
#'   `FALSE`, return table as percentages within each geography.
#'
#' @return A data frame with geographic identifier column(s) and six columns
#'   `white`, `black`, etc. containing the counts or proportion of
#'   residents in each racial group.
#'
#' @examples \dontrun{
#' census_race_geo_table("us", year=2010)
#' census_race_geo_table("state", year=2021, survey="acs1")
#' census_race_geo_table("state", year=2021, survey="acs1", GEOIDs=FALSE)
#' }
#' @concept preproc
#' @export
census_race_geo_table <- function(geo=c("us", "state", "county", "zcta", "tract"),
                                  ...,
                                  year=2010, survey=c("dec", "acs1", "acs5"),
                                  GEOIDs=TRUE, counts=TRUE) {
    rlang::check_installed("easycensus", "for downloading Census data.")
    survey = match.arg(survey)

    if (year == 2010 && survey == "dec" && geo == "zcta" && GEOIDs) {
        return(census_premade_table("GEOID", "zip_race_2010.rds", counts=counts))
    }
    if (year == 2010 && survey == "dec" && geo == "state" && GEOIDs) {
        return(census_premade_table("GEOID", "state_race_2010.rds", counts=counts))
    }

    if (survey == "dec") {
        if (!year %in% c(2010)) {
            cli_abort("Decennial census data only available for
                      {.arg year} = 2010 or 2020")
        }
        d_raw = easycensus::cens_get_dec("P5", geo, check_geo=TRUE)
    } else {
        if (year == 2020 && survey == "acs1") {
            cli_abort("No 1-year ACS data for 2020 due to the COVID-19 pandemic.")
        }
        d_raw = easycensus::cens_get_acs("B03002", geo, year=year,
                                         survey=survey, check_geo=TRUE) %>%
            filter(.data$hsplo_race_sub == "total",
                   .data$hispanic_or_latino_origin != "total") %>%
            select(-"hsplo_race_sub") %>%
            mutate(value = easycensus::get_est(.data$value))
    }

    d = d_raw %>%
        mutate(race = case_when(
                   .data$race == "total" ~ "total",
                   .data$hispanic_or_latino_origin == "hispanic or latino" ~ "hisp",
                   TRUE ~ as.character(easycensus::tidy_race(.data$race))
               ),
               race = case_when(
                   .data$race == "nhpi" ~ "asian",
                   .data$race == "two" ~ "other",
                   TRUE ~ race
               )) %>%
        group_by(.data$GEOID, .data$NAME, .data$race) %>%
        summarize(value = sum(.data$value),
                  .groups="drop") %>%
        pivot_wider_tiny(names_from="race") %>%
        # select("GEOID", "NAME", pop="total", pop_white="white",
        #        pop_black="black", pop_hisp="hisp", pop_asian="asian",
        #        pop_aian="aian", pop_other="other") %>%
        select("GEOID", "NAME", "white", "black", "hisp",
               "asian", "aian", "other", "total") %>%
        mutate(across(c(-"GEOID", -"NAME"), as.integer))

    if (isFALSE(counts)) { # normalize by population
        nc = ncol(d)
        for (i in 1:6) {
            d[, nc - i] = d[, nc - i] / d$total
        }
    }
    d$total <- NULL # delete col

    if (isTRUE(GEOIDs)) {
        select(d, -"NAME")
    } else {
        d = select(d, -"GEOID")
        colnames(d)[1] = geo
        d
    }
}

# Helper function to make an R|G table
census_premade_table = function(G_name, path, counts=FALSE) {
    d_cens = readRDS(system.file("extdata", path, package="birdie", mustWork=TRUE))
    if (!counts) { # normalize by population
        for (i in 1:6) {
            d_cens[, 2+i] = d_cens[, 2+i] / d_cens[, 2]
        }
    }
    d_cens = d_cens[, -2]
    colnames(d_cens) = c(G_name, "white", "black", "hisp", "asian", "aian", "other")
    as_tibble(d_cens)
}


# Helper function to make an R|G table
census_state_table = function(G_name, counts=FALSE) {
    d_cens = readRDS(system.file("extdata", "state_race_2010.rds",
                                 package="birdie", mustWork=TRUE))
    if (!counts) { # normalize by population
        for (i in 1:6) {
            d_cens[, 2+i] = d_cens[, 2+i] / d_cens[, 2]
        }
    }
    d_cens = d_cens[, -2]
    colnames(d_cens) = c(G_name, "white", "black", "hisp", "asian", "aian", "other")
    as_tibble(d_cens)
}


# Helper function to make an R|S table
# `counts` returns counts
# `flip` returns S|R rather than R|S
census_surname_table = function(S, S_name, counts=FALSE, flip=FALSE) {
    if (counts && flip)
        cli_abort("{.arg flip} must be {.val FALSE} if {.arg counts} is {.val TRUE}")

    d_cens = readRDS(system.file("extdata", "names_2010_counts.rds",
                                 package="birdie", mustWork=TRUE))
    if (!counts && flip) {
        p_s = rowSums(d_cens[, -1])
        p_s = p_s / sum(p_s)
        for (r in names(d_cens)[-1]) {
            d_cens[[r]] = d_cens[[r]] * p_s
            d_cens[[r]] = d_cens[[r]] / sum(d_cens[[r]])
        }
    }

    out = data.frame(last_name = unique(S)) %>%
        left_join(d_cens, by="last_name")
    missing_idx = which(is.na(out$pr_white))
    double_idx = missing_idx[is_double_name(out$last_name[missing_idx])]
    # track which indices will become <generic>
    bad = setdiff(missing_idx, double_idx)

    # double-barreled surnames
    match_1st_idx = match(stringr::word(out$last_name[double_idx], 1, 1), d_cens$last_name)
    match_2nd_idx = match(stringr::word(out$last_name[double_idx], 2, 2), d_cens$last_name)

    # neither matches
    na_ct = is.na(match_1st_idx) + is.na(match_2nd_idx)
    if (length(bad) > 0)
        bad = na.omit(c(bad, double_idx[which(na_ct == 2)]))

    # one or the other matches: use it but with half confidence
    only_1st = which(!is.na(match_1st_idx) & is.na(match_2nd_idx))
    only_2nd = which(is.na(match_1st_idx) & !is.na(match_2nd_idx))
    out[double_idx[only_1st], -1] = d_cens[match_1st_idx[only_1st], -1] / 2
    out[double_idx[only_2nd], -1] = d_cens[match_2nd_idx[only_2nd], -1] / 2

    # both match: average guesses
    both = which(na_ct == 0)
    pr_1st = as.matrix(d_cens[match_1st_idx[both], -1])
    pr_2nd = as.matrix(d_cens[match_2nd_idx[both], -1])
    sum1 = rowSums(pr_1st)
    sum2 = rowSums(pr_2nd)
    new_totals = pmin(sum1, sum2) / 2 # err on the side of less information
    out[double_idx[both], -1] =  0.5*new_totals*(pr_1st/sum1 + pr_2nd/sum2)

    # replace non-matches w/ <generic>
    if (length(bad) > 0) {
        out = rbind(out[-bad, ], tail(d_cens, 1))
    } else {
        out = rbind(out, tail(d_cens, 1))
    }
    if (!counts && !flip) {
        out[, -1] = out[, -1] / rowSums(out[, -1])
    }

    colnames(out) = c(S_name, "white", "black", "hisp", "asian", "aian", "other")
    as_tibble(out)
}

states <- data.frame(
    stringsAsFactors = FALSE,
    GEOID = c("01","02","04","05","06",
              "08","09","10","11","12","13","15","16","17","18",
              "19","20","21","22","23","24","25","26","27",
              "28","29","30","31","32","33","34","35","36","37",
              "38","39","40","41","42","44","45","46","47",
              "48","49","50","51","53","54","55","56","60","66",
              "69","72","74","78"),
    abbr = c("AL","AK","AZ","AR","CA",
             "CO","CT","DE","DC","FL","GA","HI","ID","IL","IN",
             "IA","KS","KY","LA","ME","MD","MA","MI","MN",
             "MS","MO","MT","NE","NV","NH","NJ","NM","NY","NC",
             "ND","OH","OK","OR","PA","RI","SC","SD","TN",
             "TX","UT","VT","VA","WA","WV","WI","WY","AS","GU",
             "MP","PR","UM","VI"),
    name =c("ALABAMA","ALASKA","ARIZONA",
            "ARKANSAS","CALIFORNIA","COLORADO","CONNECTICUT",
            "DELAWARE","DISTRICT OF COLUMBIA","FLORIDA","GEORGIA",
            "HAWAII","IDAHO","ILLINOIS","INDIANA","IOWA","KANSAS",
            "KENTUCKY","LOUISIANA","MAINE","MARYLAND",
            "MASSACHUSETTS","MICHIGAN","MINNESOTA","MISSISSIPPI","MISSOURI",
            "MONTANA","NEBRASKA","NEVADA","NEW HAMPSHIRE",
            "NEW JERSEY","NEW MEXICO","NEW YORK","NORTH CAROLINA",
            "NORTH DAKOTA","OHIO","OKLAHOMA","OREGON","PENNSYLVANIA",
            "RHODE ISLAND","SOUTH CAROLINA","SOUTH DAKOTA",
            "TENNESSEE","TEXAS","UTAH","VERMONT","VIRGINIA",
            "WASHINGTON","WEST VIRGINIA","WISCONSIN","WYOMING",
            "AMERICAN SAMOA","GUAM","NORTHERN MARIANA ISLANDS","PUERTO RICO",
            "U.S. MINOR OUTLYING ISLANDS","U.S. VIRGIN ISLANDS")
)
