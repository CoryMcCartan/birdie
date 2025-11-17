# National Racial Demographics

Returns the proportion of the U.S. population in six racial groups in a
given year. Group definitions necessarily follow those used by the
Census Bureau in its surname tables:

- `white`: Non-Hispanic White alone

- `black`: Non-Hispanic Black alone

- `hisp`: Hispanic, any race

- `asian`: Non-Hispanic Asian, Native Hawaiian, or Pacific Islander
  alone

- `aian`: Non-Hispanic American Indian/Alaska Native

- `other`: Non-Hispanic, two or more races, or other race

## Usage

``` r
p_r_natl(year = 2021, vap = FALSE)
```

## Arguments

- year:

  The year to return demographics for.

- vap:

  If `TRUE`, return statistics for the voting-age population (18+)
  rather than the full U.S. population.

## Value

A named numeric vector of length 6.

## Examples

``` r
p_r_natl(year=2010)
#>       white       black        hisp       asian        aian       other 
#> 0.637474968 0.122061191 0.163492546 0.048411064 0.007278155 0.021282076 
```
