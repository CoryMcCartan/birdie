#' Preprocess Last Names and Geographic Identifiers
#'
#' These functions are called automatically by [bisg()] but may be useful,
#' especially when geographic variables are included in a [birdie()] model.
#' `proc_zip()` and `proc_state()` preprocess their corresponding geographic
#' identifiers. States are partially matched to state names and abbreviations
#' and are returned as FIPS codes. ZIP codes are crosswalked to Census ZCTAs.
#' Missing identifiers are replaced with `"<none>"`.
#' `proc_name()` processes last names in accordance with Census processing rules
#' (<https://www2.census.gov/topics/genealogy/2010surnames/surnames.pdf>).
#' Names are converted to Latin characters, capitalized, stripped of prefixes
#' and suffixes, and otherwise standardized.
#'
#' @param x A character vector of names or geographic identifiers to process
#' @param to_latin If `TRUE`, convert names to Latin characters only. Strongly
#'   recommended if non-Latin characters are present, since these will not match
#'   Census tables. However, the conversion is slightly time-consuming and so
#'   can be disabled with this flag.
#'
#' @returns A processed character vector
#'
#' @examples
#' proc_name("Smith Jr.")
#' proc_zip("00501")
#' proc_state("Washington")
#' @concept preproc
#' @name preproc
NULL

#' @describeIn preproc Match ZIP codes to ZCTAs and fill in missing values.
#' @export
proc_zip = function(x) {
    coalesce(zip_xw$zcta[match(x, zip_xw$zip)], "<none>")
}

#' @describeIn preproc Match state names and abbreviations and fill in missing values.
#' @export
proc_state = function(x) {
    x = stringr::str_to_upper(x)
    idx = coalesce(
        pmatch(x, states$GEOID, duplicates.ok=TRUE),
        pmatch(x, states$abbr, duplicates.ok=TRUE),
        pmatch(x, states$name, duplicates.ok=TRUE)
    )
    coalesce(states$GEOID[idx], "<none>")
}

#' @describeIn preproc Process names to a Census-standardized format.
#' @export
proc_name = function(x, to_latin=TRUE) {
    x = stringr::str_to_upper(x)
    if (to_latin) x = stringi::stri_trans_general(x, "Latin-ASCII")
    x = if_else(str_starts(x, "[A-Z] [A-Z]$"), NA_character_, x)
    x = if_else(str_starts(x, "([A-Z] ){2,}[A-Z]$"), str_remove_all(x, " "), x)
    x = str_remove_all(x, "[.,\\'\"!@#$%^&*/?~`]")
    x = if_else(x %in% c("JUNIOR", "SENIOR", "THIRD", "CRUTHIRD"), x,
                str_remove(x, "(JUNIOR|SENIOR|THIRD|JR|III| II| IV| J R| S R)$"))
    x = str_remove(x, "^(JR|III|II |J R |S R )")
    x = if_else(stringr::str_length(x) >= 7, str_remove(x, "SR$"), x)
    x = str_replace_all(x, "^(MC|MAC|O) ", "\\1")
    x = str_replace_all(x, "-", " ")
    x = stringr::str_squish(x)
    x = if_else(stringr::str_length(x) == 2 & x != "NG" &
                    !str_detect(x, "[AEIOUY]"), NA_character_, x)
    x = if_else(x %in% c("DECLINE TO STATE", "DONT KNOW", "DECLINE",
                         "NO NAME", "NO NOMBRE", "SAME AS ABOVE"),
                NA_character_, x)
    x
}

is_double_name = function(x) {
    str_detect(x, "^\\S+ \\S+$")
}
