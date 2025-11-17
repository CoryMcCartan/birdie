# A pseudo-voterfile

A dataset containing 5,000 fake voter records. Created by randomizing a
subset of the North Carolina voter file. Turnout records are completely
randomly generated.

## Usage

``` r
pseudo_vf
```

## Format

A data frame with 5,000 rows and 4 records:

- last_name:

  Voter's last name

- zip:

  5-digit ZIP code. May be NA

- race:

  One of "white", "black", "hisp", "asian", "aian", or "other"

- turnout:

  1 if the voter voted in the most recent election, 0 otherwise

## Source

<https://www.ncsbe.gov/results-data/voter-registration-data>

## Examples

``` r
data(pseudo_vf)
print(pseudo_vf)
#> # A tibble: 5,000 × 4
#>    last_name zip   race  turnout
#>    <fct>     <fct> <fct> <fct>  
#>  1 BEAVER    28748 white yes    
#>  2 WILLIAMS  28144 black no     
#>  3 ROSEN     28270 white yes    
#>  4 SMITH     28677 black yes    
#>  5 FAY       28748 white no     
#>  6 CHURCH    28215 white yes    
#>  7 JOHNSON   28212 black yes    
#>  8 SZCZYGIEL NA    white yes    
#>  9 SUMMERS   28152 black yes    
#> 10 STARLING  28650 white yes    
#> # ℹ 4,990 more rows
```
