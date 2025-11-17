# BIRDiE Complete-Data Model Families

BIRDiE supports a number of complete-data outcome models, including
categorical regression models. Models specific to BIRDiE are listed
here. See the Details section of
[`birdie()`](http://corymccartan.com/birdie/reference/birdie.md) for
more information about each model.

## Usage

``` r
cat_dir(link = "identity")

cat_mixed(link = "softmax")
```

## Arguments

- link:

  The link function. Only one option available for categorical
  regression models.

## Value

A list of class `family` containing the specification.

## Examples

``` r
cat_dir()
#> 
#> Family: cat_dir 
#> Link function: identity 
#> 
cat_mixed()
#> 
#> Family: cat_mixed 
#> Link function: softmax 
#> 
```
