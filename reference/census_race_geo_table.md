# Download Census Race Data

Downloads and prepares race-by-geography tables from U.S. census data,
using the
[`easycensus`](http://corymccartan.com/easycensus/reference/easycensus-package.md)
package. Requires that an api key be set up through
[`easycensus::cens_auth()`](http://corymccartan.com/easycensus/reference/cens_auth.md)
in that package, usually by storing it in the `CENSUS_API_KEY`
environment variable. Supports data from the decennial census and the
American Community Survey at a variety of levels of geographic detail.
The output of this function can be used directly in
[`bisg()`](http://corymccartan.com/birdie/reference/bisg.md).

## Usage

``` r
census_race_geo_table(
  geo = c("us", "state", "county", "zcta", "tract"),
  ...,
  year = 2010,
  survey = c("dec", "acs1", "acs5"),
  GEOIDs = TRUE,
  counts = TRUE
)
```

## Arguments

- geo:

  The geographic level to return. Common options are listed in the
  function signature, but any of the geographies listed at
  [`easycensus::cens_geo()`](http://corymccartan.com/easycensus/reference/cens_geo.md)
  may be used.

- ...:

  Further subgeographies to return, as in
  [`easycensus::cens_geo()`](http://corymccartan.com/easycensus/reference/cens_geo.md).

- year:

  The year for the data

- survey:

  The data product to use: either the decennial census (`"dec"`), or the
  the 1-year or 5-year ACS.

- GEOIDs:

  If `TRUE`, return the `GEOID` column as the unique geographic
  identifier; if `FALSE`, return a human-readable name. For example,
  with `geo="state"`, setting `GEOIDs=FALSE` will return a column named
  `state` with entries like `"Massachusetts"`.

- counts:

  If `TRUE`, return the table as actual population counts; if `FALSE`,
  return table as percentages within each geography.

## Value

A data frame with geographic identifier column(s) and six columns
`white`, `black`, etc. containing the counts or proportion of residents
in each racial group.

## Examples

``` r
census_race_geo_table("zcta", year=2010)
#> # A tibble: 33,120 × 7
#>    GEOID white black  hisp asian  aian other
#>    <chr> <int> <int> <int> <int> <int> <int>
#>  1 00601    80     2 18486     1     1     0
#>  2 00602   216    13 41265    15     0    11
#>  3 00603   628   101 53877    50     2    31
#>  4 00606    32     3  6575     4     0     1
#>  5 00610   187    22 28789     8     0    10
#>  6 00612   439    45 66435    48     5    38
#>  7 00616    36    11 10965     1     1     3
#>  8 00617   110    17 24450    15     0     5
#>  9 00622   129     8  7708     3     0     5
#> 10 00623   326    28 42665    18     1    23
#> # ℹ 33,110 more rows
if (FALSE) { # \dontrun{
# Census API key required
census_race_geo_table("us", year=2010)
census_race_geo_table("state", year=2021, survey="acs1")
census_race_geo_table("county", state="NH", year=2020, GEOIDs=FALSE)
} # }
```
