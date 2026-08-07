#############
# Prepare Census name tables
library(tidyverse)
library(readxl)
library(here)
devtools::load_all(".")


# These are the complete tables linked from
# https://www.census.gov/topics/population/genealogy/data/2020_names.html.
# `edited` identifies the Census tables in which negative values from the
# disclosure-avoidance mechanism have been removed. Keep both URLs here so the
# choice between the research and edited products is explicit and reproducible.
name_table_urls = tribble(
    ~name_type, ~edited, ~url,
    "first", TRUE,
    "https://www2.census.gov/topics/genealogy/2020surnames/Names2020_FirstNames_RaceHispanic.xlsx",
    "first", FALSE,
    "https://www2.census.gov/topics/genealogy/2020surnames/Names2020_FirstNames_RaceHispanic_WithNegatives.xlsx",
    "last", TRUE,
    "https://www2.census.gov/topics/genealogy/2020surnames/Names2020_LastNames_RaceHispanic.xlsx",
    "last", FALSE,
    "https://www2.census.gov/topics/genealogy/2020surnames/Names2020_LastNames_RaceHispanic_WithNegatives.xlsx"
)

download_name_table = function(name_type=c("last", "first"), edited=TRUE,
                               dest_dir=here("data-raw/names")) {
    name_type = match.arg(name_type)
    info = filter(name_table_urls, .data$name_type == .env$name_type,
                  .data$edited == .env$edited)
    stopifnot(nrow(info) == 1L)

    dir.create(dest_dir, recursive=TRUE, showWarnings=FALSE)
    path = file.path(dest_dir, basename(info$url))
    if (!file.exists(path))
        download.file(info$url, path, mode="wb")
    path
}

# Normalize workbook labels rather than relying on column positions. The first-
# and last-name race tables have the same substantive columns, but may differ in
# title rows and sheet layout.
clean_xlsx_label = function(x) {
    x |>
        str_replace_all("%", " percent ") |>
        str_to_lower() |>
        str_replace_all("[^a-z0-9]+", " ") |>
        str_squish()
}

find_name_header = function(x) {
    score = apply(x, 1, function(row) {
        labels = clean_xlsx_label(as.character(row))
        sum(str_detect(labels, "^(name|rank|count|number)$") |
            str_detect(labels, "white|black|hispanic|asian|pacific|aian|two or more"))
    })
    which.max(score)
}

read_name_sheet = function(path, sheet) {
    raw = read_excel(path, sheet=sheet, col_names=FALSE,
                     col_types="text", .name_repair="minimal")
    header_row = find_name_header(raw)
    header = clean_xlsx_label(unlist(raw[header_row, ], use.names=FALSE))
    d = raw[-seq_len(header_row), ]
    names(d) = make.unique(header)
    select(d, where(~ any(!is.na(.x))))
}

pick_name_column = function(nms, pattern, label) {
    idx = which(str_detect(nms, regex(pattern, ignore_case=TRUE)))
    if (length(idx) != 1L)
        stop("Could not uniquely identify the ", label, " column.")
    nms[idx]
}

pick_race_column = function(nms, pattern, label, value_type) {
    idx = which(str_detect(nms, regex(pattern, ignore_case=TRUE)))
    if (value_type == "percent") {
        idx_type = idx[str_detect(nms[idx], "percent|pct|share")]
    } else {
        idx_type = idx[str_detect(nms[idx], "count|number|frequency")]
    }
    if (length(idx_type) == 1L)
        return(nms[idx_type])
    if (length(idx) == 1L)
        return(nms[idx])
    stop("Could not uniquely identify the ", label, " ", value_type, " column.")
}

# Convert a race-and-Hispanic-origin workbook to a common intermediate schema.
# Census reports mutually exclusive Hispanic and non-Hispanic race categories.
# Asian and Native Hawaiian/Other Pacific Islander are combined for consistency
# with birdie's existing surname and geography tables; multiracial is `other`.
parse_name_race_table = function(path, name_type=c("last", "first"),
                                 value_type=c("percent", "count")) {
    name_type = match.arg(name_type)
    value_type = match.arg(value_type)
    sheets = excel_sheets(path)
    parsed = map(sheets, ~ read_name_sheet(path, .x))
    d = parsed[[which.max(map_int(parsed, nrow))]]
    nms = names(d)

    asian_combined = nms[str_detect(
        nms,
        "asian.*(native hawaiian|pacific)|native hawaiian.*asian"
    )]
    if (length(asian_combined) > 1L)
        stop("Could not uniquely identify the combined Asian and NHPI column.")

    cols = c(
        name = pick_name_column(nms, "^(name|first name|last name|surname)$", "name"),
        count = pick_name_column(
            nms, "^(count|number|frequency)( count| number)?$", "count"
        ),
        white = pick_race_column(nms, "white", "White", value_type),
        black = pick_race_column(nms, "black", "Black", value_type),
        hisp = pick_race_column(nms, "^hispanic", "Hispanic", value_type),
        asian = if (length(asian_combined) == 1L) asian_combined else
            pick_race_column(nms, "asian", "Asian", value_type),
        aian = pick_race_column(nms, "american indian|alaska native|aian", "AIAN",
                                value_type),
        other = pick_race_column(nms, "two or more|multiracial", "multiracial",
                                 value_type)
    )
    if (length(asian_combined) == 0L) {
        cols = append(
            cols,
            nhpi=pick_race_column(nms, "pacific|nhpi", "NHPI", value_type),
            after=match("asian", names(cols))
        )
    }

    out = select(d, all_of(unname(cols)))
    names(out) = names(cols)
    out = out |>
        mutate(name=str_to_upper(name),
               across(-name, parse_number)) |>
        filter(!is.na(name), !is.na(count))
    if (length(asian_combined) == 1L)
        out = mutate(out, nhpi=0, .after=asian)
    structure(out, value_type=value_type)
}

# Turn reported percentages into birdie's expected integer count schema. A
# suppressed cell is allocated the unreported row total in proportion to the
# supplied population totals, following the existing 2010-table procedure.
# Percentages are rounded in the published workbook, so counts are finally
# rescaled to the published name total before integerization.
name_race_counts = function(d, p_r) {
    stopifnot(identical(names(p_r),
                        c("white", "black", "hisp", "asian", "aian", "other")))
    value_type = attr(d, "value_type")
    d = d |>
        mutate(asian = asian + nhpi) |>
        select(name, count, white, black, hisp, asian, aian, other)
    m = as.matrix(select(d, white:other))
    if (identical(value_type, "count"))
        m = m / d$count
    else
        m = m / 100

    rowtot = rowSums(m, na.rm=TRUE)
    rowwgt = as.numeric(is.na(m) %*% p_r)
    for (i in seq_len(ncol(m))) {
        fill = (p_r[i] / rowwgt) * pmax(0, 1 - rowtot)
        m[, i] = coalesce(m[, i], fill)
    }
    m = m / rowSums(m)
    m = round(m * d$count)

    # Preserve each published row total after independent rounding.
    delta = d$count - rowSums(m)
    largest = max.col(m, ties.method="first")
    m[cbind(seq_len(nrow(m)), largest)] =
        m[cbind(seq_len(nrow(m)), largest)] + delta
    storage.mode(m) = "integer"

    colnames(m) = paste0("pr_", colnames(m))
    bind_cols(tibble(name=d$name), as_tibble(m))
}

as_birdie_name_table = function(d, key) {
    generic = which(str_detect(d$name, "^ALL (OTHER|REMAINING) (FIRST |LAST )?NAMES$"))
    if (length(generic) != 1L)
        stop("Could not uniquely identify the all-other-names aggregate.")
    d$name[generic] = "<generic>"
    d = bind_rows(d[-generic, ], d[generic, ])
    rename(d, "{key}" := name)
}


# Surnames by race ------
# Use the Census-edited table in which disclosure-avoidance noise has been
# postprocessed to remove negative counts. The parser also supports the signed
# research product through `edited=FALSE` for future sensitivity analyses.
path = download_name_table("last", edited=TRUE)
d_surname = parse_name_race_table(path, "last", value_type="count")
d_surname = name_race_counts(d_surname, p_r_natl(2020)) |>
    as_birdie_name_table("last_name")
write_rds(d_surname, here("inst/extdata/names_2020_counts.rds"), compress="xz")


# First names can later use the same download, parsing, and count-conversion
# functions; only their application-facing key and output filename differ.
