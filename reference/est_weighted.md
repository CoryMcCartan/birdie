# Calculate Weighted Estimate of (Discrete) Outcomes By Race

Calculates the "standard" weighted estimator of conditional
distributions of an outcome variable \\Y\\ by race \\R\\, using BISG
probabilities. This estimator, while commonly used, is only appropriate
if \\Y \perp R \mid X, S\\, where \\S\\ and \\X\\ are the last names and
covariates (possibly including geography) used in making the BISG
probabilities. In most cases this assumption is not plausible and
[`birdie()`](http://corymccartan.com/birdie/reference/birdie.md) should
be used instead. See the references below for more discussion as to
selecting the right estimator.

Up to Monte Carlo error, the weighted estimate is equivalent to
performing multiple imputations of the race vector from the BISG
probabilities and then using them inside a weighted average or linear
regression.

## Usage

``` r
est_weighted(
  r_probs,
  formula,
  data = NULL,
  weights = NULL,
  prefix = "pr_",
  se_boot = 0
)

# S3 method for class 'est_weighted'
print(x, ...)

# S3 method for class 'est_weighted'
summary(object, ...)
```

## Arguments

- r_probs:

  A data frame or matrix of BISG probabilities, with one row per
  individual. The output of
  [`bisg()`](http://corymccartan.com/birdie/reference/bisg.md) can be
  used directly here.

- formula:

  A two-sided formula object describing the estimator structure. The
  left-hand side is the outcome variable, which must be discrete.
  Subgroups for which to calculate estimates may be specified by adding
  covariates on the right-hand side. Subgroup estimates are available
  with `coef(..., subgroup=TRUE)` and `tidy(..., subgroup=TRUE)`.

- data:

  An optional data frame containing the variables named in `formula`.

- weights:

  An optional numeric vector specifying weights.

- prefix:

  If `r_probs` is a data frame, the columns containing racial
  probabilities will be selected as those with names starting with
  `prefix`. The default will work with the output of
  [`bisg()`](http://corymccartan.com/birdie/reference/bisg.md).

- se_boot:

  The number of bootstrap replicates to use to compute an approximate
  covariance matrix for the estimator. If no bootstrapping is used, an
  analytical estimate of standard errors will be returned as `$se`. For
  bootstrapping, when there are fewer than 1,000 individuals or 100 or
  fewer replicates, a Bayesian bootstrap is used instead (i.e., weights
  are drawn from a \\\text{Dirichlet}(1, 1, ..., 1)\\ distribution,
  which produces more reliable estimates.

- ...:

  Additional arguments to generic methods (ignored).

- object, x:

  An object of class `est_weighted`.

## Value

An object of class `est_weighted`, inheriting from
[`birdie`](http://corymccartan.com/birdie/reference/birdie-class.md),
for which many methods are available. The model estimates may be
accessed with [`coef()`](https://rdrr.io/r/stats/coef.html). Uncertainty
estimates, if available, can be accessed with `$se` and
[`vcov.birdie()`](http://corymccartan.com/birdie/reference/birdie-class.md).

## Methods (by generic)

- `print(est_weighted)`: Print a summary of the model fit.

- `summary(est_weighted)`: Print a more detailed summary of the model
  fit.

## References

McCartan, C., Fisher, R., Goldin, J., Ho, D.E., & Imai, K. (2025).
Estimating Racial Disparities when Race is Not Observed. *Journal of the
American Statistical Association*. Available at
[doi:10.1080/01621459.2025.2526695](https://doi.org/10.1080/01621459.2025.2526695)
.

## Examples

``` r
data(pseudo_vf)

r_probs = bisg(~ nm(last_name) + zip(zip), data=pseudo_vf)

# Process zip codes to remove missing values
pseudo_vf$zip = proc_zip(pseudo_vf$zip)

est_weighted(r_probs, turnout ~ 1, data=pseudo_vf)
#> Weighted estimator
#> Formula: turnout ~ 1
#>    Data: pseudo_vf
#> Number of obs: 5,000; groups: 1
#> Estimated distribution:
#>     white black  hisp asian aian other
#> no  0.315 0.335 0.357 0.468 0.53 0.329
#> yes 0.685 0.665 0.643 0.532 0.47 0.671

est = est_weighted(r_probs, turnout ~ zip, data=pseudo_vf)
tidy(est, subgroup=TRUE)
#> # A tibble: 7,416 × 4
#>    zip   turnout race  estimate
#>    <chr> <chr>   <chr>    <dbl>
#>  1 28748 no      white    0.276
#>  2 28748 yes     white    0.724
#>  3 28748 no      black    0.340
#>  4 28748 yes     black    0.660
#>  5 28748 no      hisp     0.147
#>  6 28748 yes     hisp     0.853
#>  7 28748 no      asian    0.284
#>  8 28748 yes     asian    0.716
#>  9 28748 no      aian     0.211
#> 10 28748 yes     aian     0.789
#> # ℹ 7,406 more rows
```
