# Compute Racial Disparities from Model Estimates

This function lets you easily compute differences in conditional
expectations between all pairs of specified racial groups.

## Usage

``` r
disparities(x, subgroup = FALSE, races = TRUE)
```

## Arguments

- x:

  A `birdie` model object.

- subgroup:

  If `TRUE`, return subgroup-level (rather than marginal) disparity
  estimates.

- races:

  A character vector of racial groups to compute disparities for. The
  special value `TRUE`, the default, computes disparities for all racial
  groups.

## Value

A data frame containing a row with every possible disparity for the
specified `races`, which are identified by columns `race_1` and
`race_2`. The reported disparity is `estimate_1 - estimate_2`.

## Examples

``` r
data(pseudo_vf)
r_probs = bisg(~ nm(last_name) + zip(zip), data=pseudo_vf)
fit = birdie(r_probs, turnout ~ 1, data=pseudo_vf)

disparities(fit)
#> # A tibble: 60 × 6
#>    race_1 race_2 turnout estimate_1 estimate_2 disparity
#>    <chr>  <chr>  <chr>        <dbl>      <dbl>     <dbl>
#>  1 black  white  no           0.358      0.300    0.0572
#>  2 black  white  yes          0.642      0.700   -0.0572
#>  3 hisp   white  no           0.391      0.300    0.0910
#>  4 hisp   white  yes          0.609      0.700   -0.0910
#>  5 asian  white  no           0.610      0.300    0.310 
#>  6 asian  white  yes          0.390      0.700   -0.310 
#>  7 aian   white  no           0.768      0.300    0.467 
#>  8 aian   white  yes          0.232      0.700   -0.467 
#>  9 other  white  no           0.269      0.300   -0.0313
#> 10 other  white  yes          0.731      0.700    0.0313
#> # ℹ 50 more rows
disparities(fit, races=c("white", "black"))
#> # A tibble: 4 × 6
#>   race_1 race_2 turnout estimate_1 estimate_2 disparity
#>   <chr>  <chr>  <chr>        <dbl>      <dbl>     <dbl>
#> 1 black  white  no           0.358      0.300    0.0572
#> 2 black  white  yes          0.642      0.700   -0.0572
#> 3 white  black  no           0.300      0.358   -0.0572
#> 4 white  black  yes          0.700      0.642    0.0572
```
