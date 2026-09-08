## Test environments
* local R installation (macOS 14.5, arm64), R 4.5.0
* Windows Server 2022 x64 (build 20348) (on winbuilder) (devel)
* macOS 14.5 (on GitHub Actions), (release)
* Windows Server 2022 10.0.20348 (on GitHub Actions), (release)
* ubuntu 22.04.4 (on GitHub Actions), (devel)
* ubuntu 22.04.4 (on GitHub Actions), (release)
* ubuntu 22.04.4 (on GitHub Actions), (old release R 4.4)

## R CMD check results

0 errors | 0 warnings | 1 note

* This submission fixes installation errors on CRAN due to a recent update of
  StanHeaders.

* GNU make is a SystemRequirement in order to compile Stan models.

* Three examples are \dontrun{} in census_race_geo_table() since they require an 
  API key. The API key requirement is now documented in the examples, not just
  in the function documentation. This function is tested locally in
  tests/testthat/test-census.R.
  
* Examples are \donttest{} in birdie() since they generally take more than 5
  seconds to run. 
