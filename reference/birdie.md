# Fit BIRDiE Models

Fits one of three possible Bayesian Instrumental Regression for
Disparity Estimation (BIRDiE) models to BISG probabilities and
covariates. The simplest Categorical-Dirichlet model
([`cat_dir()`](http://corymccartan.com/birdie/reference/birdie-family.md))
is appropriate when there are no covariates or when all covariates are
discrete and fully interacted with another. The more general Categorical
mixed-effects model
([`cat_mixed()`](http://corymccartan.com/birdie/reference/birdie-family.md))
is a supports any number of fixed effects and up to one random
intercept. For continuous outcomes a Normal linear model is available
([`gaussian()`](https://rdrr.io/r/stats/family.html)).

## Usage

``` r
birdie(
  r_probs,
  formula,
  data,
  family = cat_dir(),
  prior = NULL,
  weights = NULL,
  algorithm = c("em", "gibbs", "em_boot"),
  iter = 400,
  warmup = 50,
  prefix = "pr_",
  ctrl = birdie.ctrl()
)
```

## Arguments

- r_probs:

  A data frame or matrix of BISG probabilities, with one row per
  individual. The output of
  [`bisg()`](http://corymccartan.com/birdie/reference/bisg.md) can be
  used directly here.

- formula:

  A two-sided formula object describing the model structure. The
  left-hand side is the outcome variable, which must be discrete. A
  single random intercept term, denoted with a vertical bar
  (`"(1 | <term>)"`), is supported on the right-hand side.

- data:

  An optional data frame containing the variables named in `formula`.

- family:

  A description of the complete-data model type to fit. Options are:

  - [`cat_dir()`](http://corymccartan.com/birdie/reference/birdie-family.md):
    Categorical-Dirichlet model. All covariates must be fully
    interacted.

  - [`cat_mixed()`](http://corymccartan.com/birdie/reference/birdie-family.md):
    Categorical mixed-effects model. Up to one random effect is
    supported.

  - [`gaussian()`](https://rdrr.io/r/stats/family.html): Linear model.

  See the Details section below for more information on the various
  models.

- prior:

  A list with entries specifying the model prior.

  - For the `cat_dir` model, the only entry is `alpha`, which should be
    a matrix of Dirichlet hyperparameters. The matrix should have one
    row for every level of the outcome variable and one column for every
    racial group. The default prior (used when `prior=NULL`) is an
    empirical Bayes prior equal to the weighted-mean estimate of the
    outcome-race table. A fully noninformative prior with all entries
    set to \\\epsilon\\ can be obtained by setting `prior=NA`. When
    `prior=NULL` and `algorithm="em"` or `"em_boot"`, 1 is added to the
    prior so that the posterior mode, rather than the mean, is shrunk
    toward these values.

  - For the `cat_mixed` model, the `prior` list should contain three
    scalar entries: `scale_int`, the standard deviation on the Normal
    prior for the intercepts (which control the global estimates of
    `Y|R`), `scale_beta`, the standard deviation on the Normal prior for
    the fixed effects, and `scale_sigma`, the prior mean of the standard
    deviation of the random intercepts. These can be a single scalar or
    a vector with an entry for each racial group.

  - For the `gaussian` model, the `prior` list should contain two
    entries: `scale_int`, controlling, the standard deviation on the
    Normal prior for the intercepts (which control the global estimates
    of `Y|R`), and `scale_beta`, controlling the standard deviation on
    the Normal prior for the fixed effects. These must be a single
    scalar. Each is expressed in terms of the estimated residual
    standard deviation (i.e., they are multiplied together to form the
    "true" prior).

  The prior is stored after model fitting in the `$prior` element of the
  fitted model object.

- weights:

  An optional numeric vector specifying likelihood weights.

- algorithm:

  The inference algorithm to use. One of 3 options:

  - `"em"`: An expectation-maximization algorithm which will perform
    inference for the maximum a posteriori (MAP) parameter values.
    Computationally efficient and supported by all the model families.
    No uncertainty quantification.

  - `"gibbs"`: A Gibbs sampler for performing full Bayesian inference.
    Generally more computationally demanding than the EM algorithm, but
    provides uncertainty quantification. Currently supported for
    [`cat_dir()`](http://corymccartan.com/birdie/reference/birdie-family.md)
    and [`gaussian()`](https://rdrr.io/r/stats/family.html) model
    families. Computation-reliability tradeoff can be controlled with
    `iter` argument.

  - `"em_boot"`: Bootstrapped version of EM algorithm. Number of
    bootstrap replicates controlled by `iter` parameter. Provides
    approximate uncertainty quantification. Currently supported for
    [`cat_dir()`](http://corymccartan.com/birdie/reference/birdie-family.md)
    and [`gaussian()`](https://rdrr.io/r/stats/family.html) model
    families.

- iter:

  The number of post-warmup Gibbs samples, or the number of bootstrap
  replicates to use to compute approximate standard errors for the main
  model estimates. Only available when `family=cat_dir()` or
  [`gaussian()`](https://rdrr.io/r/stats/family.html). Ignored if
  `algorithm="em"`.

  For bootstrapping, when there are fewer than 1,000 individuals or 100
  or fewer replicates, a Bayesian bootstrap is used instead (i.e.,
  weights are drawn from a \\\text{Dirichlet}(1, 1, ..., 1)\\
  distribution, which produces more reliable estimates.

- warmup:

  Number of warmup iterations for Gibbs sampling. Ignored unless
  `algorithm="gibbs"`.

- prefix:

  If `r_probs` is a data frame, the columns containing racial
  probabilities will be selected as those with names starting with
  `prefix`. The default will work with the output of
  [`bisg()`](http://corymccartan.com/birdie/reference/bisg.md).

- ctrl:

  A list containing control parameters for the EM algorithm and
  optimization routines. A list in the proper format can be made using
  [`birdie.ctrl()`](http://corymccartan.com/birdie/reference/birdie.ctrl.md).

## Value

An object of class
[`birdie`](http://corymccartan.com/birdie/reference/birdie-class.md),
for which many methods are available. The model estimates may be
accessed with
[`coef.birdie()`](http://corymccartan.com/birdie/reference/birdie-class.md),
and updated BISG probabilities (conditioning on the outcome) may be
accessed with
[`fitted.birdie()`](http://corymccartan.com/birdie/reference/birdie-class.md).
Uncertainty estimates, if available, can be accessed with `$se` and
[`vcov.birdie()`](http://corymccartan.com/birdie/reference/birdie-class.md).

## Details

By default, `birdie()` uses an expectation-maximization (EM) routine to
find the maximum *a posteriori* (MAP) estimate for the specified model.
Asymptotic variance-covariance matrices for the MAP estimate are
available for the Categorical-Dirichlet and Normal linear models via
bootstrapping. Full Bayesian inference is supported via Gibbs sampling
for the Categorical-Dirichlet and Normal linear models as well.

Whatever model or method is used, a finite-population estimate of the
outcome-given-race distribution for the entire observed sample is always
calculated and stored as `$est` in the returned object, which can be
accessed with
[`coef.birdie()`](http://corymccartan.com/birdie/reference/birdie-class.md)
as well.

The Categorical-Dirichlet model is specified as follows: \$\$ Y_i \mid
R_i, X_i, \Theta \sim \text{Categorical}(\theta\_{R_iX_i}) \\
\theta\_{rx} \sim \text{Dirichlet}(\alpha_r), \$\$ where \\Y\\ is the
outcome variable, \\R\\ is race, \\X\\ are covariates (fixed effects),
and \\\theta\_{rx}\\ and \\\alpha_r\\ are vectors with length matching
the number of levels of the outcome variable. There is one vector
\\\theta\_{rx}\\ for every combination of race and covariates, hence the
need for `formula` to either have no covariates or a fully interacted
structure.

The Categorical mixed-effects model is specified as follows: \$\$ Y_i
\mid R_i, X_i, \Theta \sim \text{Categorical}(g^{-1}(\mu\_{R_iX_i})) \\
\mu\_{rxy} = W\beta\_{ry} + Zu\_{ry} \\ u\_{r} \mid \vec\sigma\_{r}, L_r
\sim \mathcal{N}(0,
\text{diag}(\vec\sigma\_{r})C_r\text{diag}(\vec\sigma\_{r})) \\
\beta\_{ry} \sim \mathcal{N}(0, s^2\_{r\beta}) \\ \sigma\_{ry} \sim
\text{Inv-Gamma}(4, 3s\_{r\sigma}) \\ C_r \sim \text{LKJ}(2), \$\$ where
\\\beta\_{ry}\\ are the fixed effects, \\u\_{ry}\\ is the random
intercept, and \\g\\ is a softmax link function. Estimates for
\\\beta\_{ry}\\ and \\\sigma\_{ry}\\ are stored in the `$beta` and
`$sigma` elements of the fitted model object.

The Normal linear model is specified as follows: \$\$ Y_i \mid R_i, \vec
X_i, \Theta \sim \mathcal{N}(\vec X_i^\top\vec\theta, \sigma^2) \\
\sigma^2 \sim \text{Inv-Gamma}(n\_\sigma/2, l\_\sigma^2 n\_\sigma/2) \\
\beta\_{\text{intercept}} \sim \mathcal{N}(0, s^2\_\text{int}) \\
\beta_k \sim \mathcal{N}(0, s^2\_\beta), \\ \$\$ where \\\vec\theta\\ is
a vector of linear model coefficients. Estimates for \\\theta\\ and
\\\sigma\\ are stored in the `$beta` and `$sigma` elements of the fitted
model object.

More details on the models and their properties may be found in the
paper referenced below.

## References

McCartan, C., Fisher, R., Goldin, J., Ho, D.E., & Imai, K. (2025).
Estimating Racial Disparities when Race is Not Observed. *Journal of the
American Statistical Association*. Available at
[doi:10.1080/01621459.2025.2526695](https://doi.org/10.1080/01621459.2025.2526695)
.

## Examples

``` r
# \donttest{
data(pseudo_vf)

r_probs = bisg(~ nm(last_name) + zip(zip), data=pseudo_vf)

# Process zip codes to remove missing values
pseudo_vf$zip = proc_zip(pseudo_vf$zip)

fit = birdie(r_probs, turnout ~ 1, data=pseudo_vf)
#> Using weakly informative empirical Bayes prior for Pr(Y | R)
#> This message is displayed once every 8 hours.
print(fit)
#> Categorical-Dirichlet BIRDiE model
#> Formula: turnout ~ 1
#>    Data: pseudo_vf
#> Number of obs: 5,000
#> Estimated distribution:
#>     white black  hisp asian  aian other
#> no    0.3 0.358 0.391  0.61 0.768 0.269
#> yes   0.7 0.642 0.609  0.39 0.232 0.731
fit$se # uncertainty quantification
#> NULL

fit = birdie(r_probs, turnout ~ zip, data=pseudo_vf, algorithm="gibbs")

fit = birdie(r_probs, turnout ~ (1 | zip), data=pseudo_vf,
             family=cat_mixed(), ctrl=birdie.ctrl(abstol=1e-3))
#> Using default prior for Pr(Y | R):
#> → Prior scale on intercepts: 2.0
#> → Prior scale on fixed effects coefficients: 0.2
#> → Prior mean of random effects standard deviation: 0.10
#> This message is displayed once every 8 hours.
#> ⠙ EM iterations 3 done (1.2/s) | 2.6s
#> ⠹ EM iterations 5 done (0.9/s) | 5.5s
#> ⠸ EM iterations 8 done (1/s) | 7.7s
#> ⠸ EM iterations 9 done (0.96/s) | 9.4s

summary(fit)
#> Categorical mixed-effects BIRDiE model
#> Formula: turnout ~ (1 | zip)
#>    Data: pseudo_vf
#> 
#> 9 iterations and 9.4 secs to convergence
#> 
#> Number of observations: 5,000
#> Number of groups: 618
#> 
#> Entropy decrease from marginal race distribution:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> -0.4795  0.1505  0.4227  0.4472  0.7798  1.0721 
#> Entropy decrease from BISG probabilities:
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> -0.85529 -0.02126  0.01974  0.06492  0.13638  1.16925 
#> 
#> Estimated outcome-by-race distribution:
#>     white black  hisp asian aian other
#> no  0.296 0.366 0.376 0.559 0.63 0.356
#> yes 0.704 0.634 0.624 0.441 0.37 0.644
coef(fit)
#>         white     black      hisp     asian      aian     other
#> no  0.2955088 0.3655487 0.3756242 0.5590668 0.6296485 0.3561447
#> yes 0.7044912 0.6344513 0.6243758 0.4409332 0.3703515 0.6438553
fitted(fit)
#> # A tibble: 5,000 × 6
#>    pr_white pr_black   pr_hisp   pr_asian  pr_aian pr_other
#>       <dbl>    <dbl>     <dbl>      <dbl>    <dbl>    <dbl>
#>  1  0.952   0.000953 0.0141    0.000603   0.0111     0.0213
#>  2  0.00190 0.987    0.000289  0.000162   0.000559   0.0105
#>  3  0.965   0.00669  0.000142  0.00000522 0.000423   0.0281
#>  4  0.596   0.368    0.00194   0.0000264  0.00124    0.0325
#>  5  0.983   0.00352  0.0000961 0.00191    0.000397   0.0106
#>  6  0.555   0.271    0.0985    0.00637    0.00269    0.0658
#>  7  0.117   0.773    0.0480    0.00545    0.00169    0.0548
#>  8  0.964   0.00421  0.0313    0          0          0     
#>  9  0.739   0.217    0.00995   0.000770   0.000795   0.0322
#> 10  0.864   0.0979   0.0139    0.000968   0.00145    0.0214
#> # ℹ 4,990 more rows
# }
```
