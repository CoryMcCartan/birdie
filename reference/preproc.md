# Preprocess Last Names and Geographic Identifiers

These functions are called automatically by
[`bisg()`](http://corymccartan.com/birdie/reference/bisg.md) but may be
useful, especially when geographic variables are included in a
[`birdie()`](http://corymccartan.com/birdie/reference/birdie.md) model.
`proc_zip()` and `proc_state()` preprocess their corresponding
geographic identifiers. States are partially matched to state names and
abbreviations and are returned as FIPS codes. ZIP codes are crosswalked
to Census ZCTAs. Missing identifiers are replaced with `"<none>"`.
`proc_name()` processes last names in accordance with Census processing
rules
(<https://www2.census.gov/topics/genealogy/2010surnames/surnames.pdf>).
Names are converted to Latin characters, capitalized, stripped of
prefixes and suffixes, and otherwise standardized.

## Usage

``` r
proc_zip(x)

proc_state(x)

proc_name(x, to_latin = TRUE)
```

## Arguments

- x:

  A character vector of names or geographic identifiers to process

- to_latin:

  If `TRUE`, convert names to Latin characters only. Strongly
  recommended if non-Latin characters are present, since these will not
  match Census tables. However, the conversion is slightly
  time-consuming and so can be disabled with this flag.

## Value

A processed character vector

## Functions

- `proc_zip()`: Match ZIP codes to ZCTAs and fill in missing values.

- `proc_state()`: Match state names and abbreviations and fill in
  missing values.

- `proc_name()`: Process names to a Census-standardized format.

## Examples

``` r
proc_name("Smith Jr.")
#> [1] "SMITH"
proc_zip("00501")
#> [1] "11742"
proc_state("Washington")
#> [1] "53"
```
