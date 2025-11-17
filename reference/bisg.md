# Bayesian Improved Surname Geocoding (BISG)

Calculates individual probabilities of belonging to racial groups given
last name, location, and other covariates (optional). The standard
function `bisg()` treats the input tables as fixed. An alternative
function `bisg_me()`, assumes that the input tables are subject to
measurement error, and uses a Gibbs sampler to impute the individual
race probabilities, using the model of Imai et al. (2022).

## Usage

``` r
bisg(
  formula,
  data = NULL,
  p_r = p_r_natl(),
  p_rgx = NULL,
  p_rs = NULL,
  save_rgx = TRUE,
  save_sr = FALSE
)

bisg_me(
  formula,
  data = NULL,
  p_r = p_r_natl(),
  p_rgx = NULL,
  p_rs = NULL,
  iter = 1000,
  warmup = 100,
  cores = 1L
)

# S3 method for class 'bisg'
summary(object, p_r = NULL, ...)

# S3 method for class 'bisg'
predict(object, adj = NULL, ...)

# S3 method for class 'bisg'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- formula:

  A formula specifying the BISG model. Must include the special term
  `nm()` to identify the surname variable. Certain geographic variables
  can be identified similarly:
  [`zip()`](https://rdrr.io/r/utils/zip.html) for ZIP codes, and
  `state()` for states. If no other predictor variables are provided,
  then `bisg()` will automatically be able to build a table of census
  data to use in inference. If other predictor variables are included,
  or if other geographic identifiers are used, then the user must
  specify the `p_rgx` argument below. The left-hand side of the formula
  is ignored. See the examples section below for sample formulas.

- data:

  The data frame containing the variables in `formula`.

- p_r:

  The prior distribution of race in the sample, as a numeric vector.
  Defaults to U.S. demographics as provided by
  [`p_r_natl()`](http://corymccartan.com/birdie/reference/p_r_natl.md).
  Can also set `p_r="est"` or `"estimate"` to estimate this from the
  geographic distribution. Since the prior distribution on race strongly
  affects the calibration of the BISG probabilities and thus the
  accuracy of downstream estimates, users are encouraged to think
  carefully about an appropriate value for `p_r`. If no prior
  information on the racial makeup of the sample is available, and yet
  the sample is very different from the overall U.S. population, then
  `p_r="estimate"` will likely produce superior results.

- p_rgx:

  The distribution of race given location (G) and other covariates (X)
  specified in `formula`. Should be provided as a data frame, with
  columns matching the predictors in `formula`, and additional columns
  for each racial group containing the conditional probability for that
  racial group given the predictors. For example, if Census tracts are
  the only predictors, `p_rgx` should be a data frame with a tract
  column and columns `white`, `black`, etc. containing the racial
  distribution of each tract. If `formula` contains only labeled terms
  (like [`zip()`](https://rdrr.io/r/utils/zip.html)), then by default
  `p_rgx` will be constructed automatically from the most recent Census
  data. This table will be normalized by row, so it can be provided as
  population counts as well. Counts are required for `bisg_me()`. The
  [`census_race_geo_table()`](http://corymccartan.com/birdie/reference/census_race_geo_table.md)
  function can be helpful to prepare tables, as can be the `build_dec()`
  and `build_acs()` functions in the `censable` package.

- p_rs:

  The distribution of race given last name. As with `p_rgx`, should be
  provided as a data frame, with a column of names and additional
  columns for each racial group. Users should not have to specify this
  argument in most cases, as the table will be built from published
  Census surname tables automatically. Counts are required for
  `bisg_me()`. One of the last names can be `"<generic>"`, which, if
  included, will be used for any names not found in the table.

- save_rgx:

  If `TRUE`, save the `p_rgx` table (matched to each individual) as the
  `"p_rgx"` and `"gx"` attributes of the output. Necessary for some
  sensitivity analyses.

- save_sr:

  If `TRUE`, save the `p_sr` table (surname given race; matched to each
  individual as the `"p_sr"` and `"s"` attributes of the output.

- iter:

  How many sampling iterations in the Gibbs sampler

- warmup:

  How many burn-in iterations in the Gibbs sampler

- cores:

  How many parallel cores to use in computation. Around 4 seems to be
  optimal, even if more are available.

- object:

  An object of class `bisg`, the result of running `bisg()`.

- ...:

  Additional arguments to generic methods (ignored).

- adj:

  A point in the simplex that describes how BISG probabilities will be
  thresholded to produce point predictions. The probabilities are
  divided by `adj`, then the racial category with the highest
  probability is predicted. Can be used to trade off types of prediction
  error. Must be nonnegative but will be normalized to sum to 1. The
  default is to make no adjustment.

- nsim:

  The number of vectors to simulate. Defaults to 1.

- seed:

  Used to seed the random number generator. See
  [`stats::simulate()`](https://rdrr.io/r/stats/simulate.html).

## Value

An object of class `bisg`, which is just a data frame with some
additional attributes. The data frame has rows matching the input data
and columns for the race probabilities.

## Methods (by generic)

- `summary(bisg)`: Summarize predicted race probabilities. Returns
  vector of individual entropies.

- `predict(bisg)`: Create point predictions of individual race. Returns
  factor vector of individual race labels. Strongly not recommended for
  any kind of inferential purpose, as biases may be extreme and in
  unpredictable directions.

- `simulate(bisg)`: Simulate race from the `Pr(R | G, X, S)`
  distribution.

## Functions

- `bisg()`: The standard BISG model.

- `bisg_me()`: The measurement error BISG model.

## Surname Matching

The Census surname table can be inspected with the following code:

    readRDS(system.file("extdata", "names_2010_counts.rds", package="birdie"))

Surnames are processed with
[`proc_name()`](http://corymccartan.com/birdie/reference/preproc.md)
before being matched to the table. The code also recognizes
double-barrelled (hyphenated) surnames and attempts to match each part
if the overall name is not found in the surname table. Specifying
`save_rs=TRUE` will save the matched surname table and a lookup vector
that matches each individual to their surname table row. The overall
match rate is reported as part of the
[`summary()`](https://rdrr.io/r/base/summary.html) output.

## References

Elliott, M. N., Fremont, A., Morrison, P. A., Pantoja, P., and Lurie, N.
(2008). A new method for estimating race/ethnicity and associated
disparities where administrative records lack self-reported
race/ethnicity. *Health Services Research*, 43(5p1):1722–1736.

Fiscella, K. and Fremont, A. M. (2006). Use of geocoding and surname
analysis to estimate race and ethnicity. *Health Services Research*,
41(4p1):1482–1500.

Imai, K., Olivella, S., & Rosenman, E. T. (2022). Addressing census data
problems in race imputation via fully Bayesian Improved Surname
Geocoding and name supplements. *Science Advances*, 8(49), eadc9824.

## Examples

``` r
data(pseudo_vf)
bisg(~ nm(last_name), data=pseudo_vf)
#> # A tibble: 5,000 × 6
#>    pr_white pr_black pr_hisp pr_asian pr_aian pr_other
#>       <dbl>    <dbl>   <dbl>    <dbl>   <dbl>    <dbl>
#>  1    0.826  0.0701   0.0305  0.00353 0.0208    0.0491
#>  2    0.362  0.547    0.0298  0.00315 0.00454   0.0533
#>  3    0.918  0.00887  0.0397  0.00697 0.00223   0.0240
#>  4    0.620  0.293    0.0317  0.00379 0.00545   0.0460
#>  5    0.892  0.0237   0.0413  0.00831 0.00295   0.0322
#>  6    0.844  0.0790   0.0309  0.00384 0.00542   0.0370
#>  7    0.491  0.419    0.0297  0.00390 0.00548   0.0512
#>  8    0.982  0.00574  0.0120  0       0         0     
#>  9    0.713  0.194    0.0319  0.00422 0.00695   0.0497
#> 10    0.593  0.337    0.0262  0.00278 0.00406   0.0368
#> # ℹ 4,990 more rows

r_probs = bisg(~ nm(last_name) + zip(zip), data=pseudo_vf)
summary(r_probs)
#> BISG individual race probabilities
#> 
#> Unmatched surnames: 246 (4.92%)
#> Unmatched geographies/covariates: 687 (13.7%)
#> 
#> Implied marginal race distribution:
#> pr_white pr_black  pr_hisp pr_asian  pr_aian pr_other 
#>    0.641    0.215    0.074    0.020    0.007    0.043 
#> 
#> Entropy decrease from marginal distribution:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> -0.5137  0.1193  0.3315  0.3665  0.6386  1.0563 
head(predict(r_probs))
#> [1] white black white white white white
#> Levels: white black hisp asian aian other

data(pseudo_vf)
bisg_me(~ nm(last_name) + zip(zip), data=pseudo_vf)
#> # A tibble: 5,000 × 6
#>    pr_white pr_black pr_hisp pr_asian pr_aian pr_other
#>       <dbl>    <dbl>   <dbl>    <dbl>   <dbl>    <dbl>
#>  1    0.963   0.005   0.004     0.001   0.016   0.011 
#>  2    0.205   0.772   0.009     0       0.001   0.013 
#>  3    0.977   0.001   0.01      0.006   0       0.006 
#>  4    0.649   0.314   0.015     0       0.002   0.0200
#>  5    0.994   0.001   0.002     0.001   0.001   0.001 
#>  6    0.640   0.281   0.0370    0.011   0.005   0.0260
#>  7    0.180   0.747   0.0300    0.005   0.001   0.0370
#>  8    0.985   0.002   0.013     0       0       0     
#>  9    0.799   0.175   0.006     0       0.004   0.016 
#> 10    0.911   0.0750  0.005     0       0.002   0.007 
#> # ℹ 4,990 more rows
```
