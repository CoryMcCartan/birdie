
<!-- README.md is generated from README.Rmd. Please edit that file -->

# **BIRDiE**: Estimating disparities when race is not observed <img src="man/figures/logo.png" align="right" height="156" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/CoryMcCartan/birdie/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/CoryMcCartan/birdie/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Bayesian Instrumental Regression for Disparity Estimation (BIRDiE) is a
class of Bayesian models for accurately estimating conditional
distributions by race, using Bayesian Improved Surname Geocoding (BISG)
probability estimates of individual race. This package implements BIRDiE
as described in [McCartan, Goldin, Ho and Imai
(2022)](https://arxiv.org/abs/2303.02580). It also implements standard
BISG and an improved measurement-error BISG model as described in [Imai,
Olivella, and Rosenman
(2022)](https://www.science.org/doi/full/10.1126/sciadv.adc9824).

<img src="man/figures/poster.svg" style="width: 100%" alt="BIRDiE Overview Poster" />

## Installation

You can install the development version of birdie from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("CoryMcCartan/birdie")
```

## Basic Usage

A basic analysis has two steps. First, you compute BISG probability
estimates with the `bisg()` or `bisg_me()` functions. Then, you estimate
the distribution of an outcome variable by race using the `birdie()`
function.

``` r
library(birdie)

data(pseudo_vf)

head(pseudo_vf)
#> # A tibble: 6 × 4
#>   last_name zip   race  turnout
#>   <fct>     <fct> <fct> <fct>  
#> 1 BEAVER    28748 white yes    
#> 2 WILLIAMS  28144 black no     
#> 3 ROSEN     28270 white yes    
#> 4 SMITH     28677 black yes    
#> 5 FAY       28748 white no     
#> 6 CHURCH    28215 white yes
```

To compute BISG probabilities, you provide the last name and
(optionally) geography variables as part of a formula.

``` r
r_probs = bisg(~ nm(last_name) + zip(zip), data=pseudo_vf)

head(r_probs)
#> # A tibble: 6 × 6
#>   pr_white pr_black pr_hisp pr_asian  pr_aian pr_other
#>      <dbl>    <dbl>   <dbl>    <dbl>    <dbl>    <dbl>
#> 1    0.956  0.00371  0.0103 0.000674 0.00886    0.0202
#> 2    0.162  0.795    0.0122 0.00102  0.000873   0.0292
#> 3    0.943  0.00378  0.0218 0.0107   0.000386   0.0202
#> 4    0.569  0.365    0.0302 0.00114  0.00108    0.0339
#> 5    0.971  0.00118  0.0131 0.00149  0.00118    0.0125
#> 6    0.524  0.315    0.0909 0.00598  0.00255    0.0610
```

Computing regression estimates requires specifying a model structure.
Here, we’ll let the relationship between turnout and race vary by ZIP
code. This is the “no-pooling” model from McCartan et al.

``` r
fit = birdie(r_probs, turnout ~ proc_zip(zip), data=pseudo_vf)
#> Using c(1+ε, 1+ε, ..., 1+ε) prior for Pr(X | R)
#> This message is displayed once per session.

print(fit)
#> Multinomial-Dirichlet BIRDiE model
#> Formula: turnout ~ proc_zip(zip)
#>    Data: pseudo_vf
#> Number of obs: 5,000; groups: 618
#> Estimated distribution:
#>     white black  hisp asian  aian other
#> no  0.286 0.358 0.376  0.55 0.644 0.534
#> yes 0.714 0.642 0.624  0.45 0.356 0.466
```

The `proc_zip()` function fills in missing ZIP codes, among other
things. We can extract the estimated conditional distributions with
`coef()`. We can also get updated BISG probabilities that additionally
condition on turnout using `fitted()`. Additional functions allow us to
extract a tidy version of our estimates (`tidy()`) and visualize the
estimated distributiosn (`plot()`).

``` r
coef(fit)
#>         white     black      hisp     asian      aian     other
#> no  0.2855931 0.3583592 0.3761561 0.5501976 0.6443233 0.5341151
#> yes 0.7144069 0.6416408 0.6238439 0.4498024 0.3556767 0.4658849

head(fitted(fit))
#> # A tibble: 6 × 6
#>   pr_white pr_black  pr_hisp pr_asian  pr_aian pr_other
#>      <dbl>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl>
#> 1 9.46e- 1 3.68e-15 1.41e- 2 3.81e-14 1.21e- 2 2.77e- 2
#> 2 2.77e-15 9.99e- 1 5.45e-15 5.35e-15 1.25e- 3 1.99e-14
#> 3 9.52e- 1 7.52e- 3 7.34e-15 4.13e-16 1.95e-14 4.01e- 2
#> 4 6.18e- 1 3.81e- 1 3.54e-15 3.26e-13 1.56e- 3 3.02e-13
#> 5 9.90e- 1 4.33e- 3 3.00e-11 5.48e- 3 1.10e-14 1.83e-13
#> 6 5.63e- 1 2.65e- 1 9.76e- 2 6.42e- 3 2.74e- 3 6.55e- 2

tidy(fit)
#> # A tibble: 12 × 3
#>    turnout race  estimate
#>    <chr>   <chr>    <dbl>
#>  1 no      white    0.286
#>  2 yes     white    0.714
#>  3 no      black    0.358
#>  4 yes     black    0.642
#>  5 no      hisp     0.376
#>  6 yes     hisp     0.624
#>  7 no      asian    0.550
#>  8 yes     asian    0.450
#>  9 no      aian     0.644
#> 10 yes     aian     0.356
#> 11 no      other    0.534
#> 12 yes     other    0.466

plot(fit)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

A more detailed introduction to the method and software package can be
found on the [Get
Started](https://corymccartan.com/birdie/articles/birdie.html) page.
