# Class "birdie" of BIRDiE Models

The output of
[`birdie()`](http://corymccartan.com/birdie/reference/birdie.md) is an
object of class `birdie`, which supports many generic functions. Notably
`coef.birdie()` returns the main model estimates of outcome given race,
and `fitted.birdie()` returns a table analogous to the output of
[`bisg()`](http://corymccartan.com/birdie/reference/bisg.md) with
updated race probabilities.

## Usage

``` r
# S3 method for class 'birdie'
coef(object, subgroup = FALSE, ...)

# S3 method for class 'birdie'
fitted(object, ...)

# S3 method for class 'birdie'
residuals(object, x_only = FALSE, ...)

# S3 method for class 'birdie'
predict(object, adj = NULL, ...)

# S3 method for class 'birdie'
simulate(object, nsim = 1, seed = NULL, ...)

# S3 method for class 'birdie'
plot(x, log = FALSE, ...)

# S3 method for class 'birdie'
tidy(x, subgroup = FALSE, ...)

# S3 method for class 'birdie'
glance(x, ...)

# S3 method for class 'birdie'
augment(x, data, ...)

# S3 method for class 'birdie'
formula(x, ...)

# S3 method for class 'birdie'
family(object, ...)

# S3 method for class 'birdie'
nobs(object, ...)

# S3 method for class 'birdie'
vcov(object, ...)

# S3 method for class 'birdie'
print(x, ...)

# S3 method for class 'birdie'
summary(object, ...)
```

## Arguments

- object, x:

  A `birdie` model object

- subgroup:

  If `TRUE`, return subgroup-level (rather than marginal) coefficient
  estimates as a 3D array.

- ...:

  Potentially further arguments passed from other methods

- x_only:

  if `TRUE`, calculate fitted values using covariates only (i.e.,
  without using surnames).

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

- log:

  If `TRUE`, plot estimated probabilities on a log scale.

- data:

  A data frame to augment with `Pr(R | Y, X, S)` probabilities

## Value

Varies, depending on the method. See generic functions' documentation
for details.

## Details

The internal structure of `birdie` objects is not designed to be
accessed directly. The generics listed here should be used instead.

## Functions

- `coef(birdie)`: Return estimated outcome-given-race distributions.
  When `subgroup=FALSE` this always returns a finite-population estimate
  of the outcome-given-race distribution for the observed sample.

- `fitted(birdie)`: Return an updated race probability table.
  [`bisg()`](http://corymccartan.com/birdie/reference/bisg.md) estimates
  `Pr(R | G, X, S)`; this table is `Pr(R | Y, G, X, S, Theta-hat)`.

- `residuals(birdie)`: Return the residuals for the outcome variable as
  a matrix. Useful in sensitivity analyses and to get an idea of how
  well race, location, names, etc. predict the outcome.

- `predict(birdie)`: Create point predictions of individual race.
  Returns factor vector of individual race labels. Strongly not
  recommended for any kind of inferential purpose, as biases may be
  extreme and in unpredictable directions.

- `simulate(birdie)`: Simulate race from the posterior distribution
  `Pr(R | Y, G, X, S, Theta-hat)`. Does not account for uncertainty in
  model parameters.

- `plot(birdie)`: Visualize the estimated conditional distributions for
  a BIRDiE model. If available, marginal standard error estimates
  (`$se`) will be visualized with 95% confidence-level error bars.

- `tidy(birdie)`: Put BIRDiE model coefficients in a tidy format.

- `glance(birdie)`: Glance at a BIRDiE model.

- `augment(birdie)`: Augment data with individual race predictions from
  a BIRDiE model.

- `formula(birdie)`: Extract the formula used to specify a BIRDiE model.

- `family(birdie)`: Return the BIRDiE complete-data model family.

- `nobs(birdie)`: Return the number of observations used to fit a BIRDiE
  model.

- `vcov(birdie)`: Return the estimated variance-covariance matrix for
  the BIRDiE model estimates, if available.

- `print(birdie)`: Print a summary of the model fit.

- `summary(birdie)`: Print a more detailed summary of the model fit.

## Examples

``` r
methods(class="birdie")
#>  [1] augment   coef      family    fitted    formula   glance    nobs     
#>  [8] plot      predict   print     residuals simulate  summary   tidy     
#> [15] vcov     
#> see '?methods' for accessing help and source code
```
