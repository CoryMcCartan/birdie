## Test environments
* local R installation (macOS 14.5, arm64), R 4.4.0
* macOS 14.5 (on GitHub Actions), (release)
* Windows Server 2022 10.0.20348 (on GitHub Actions), (release)
* ubuntu 22.04.4 (on GitHub Actions), (devel)
* ubuntu 22.04.4 (on GitHub Actions), (release)
* ubuntu 22.04.4 (on GitHub Actions), (old release R 4.3.3)

## R CMD check results

0 errors | 0 warnings | 2 notes

* This is a new release.

* GNU make is a SystemRequirement in order to compile Stan models.

* C++17 is a requirement in order to compile Stan models as well.
  This is noted in src/Makevars and in SystemRequirements.

* Examples are \dontrun{} in census_race_geo_table() since they require an API
  key and may take some time to run. This function is tested locally in 
  tests/testthat/test-census.R.
* Examples are \donttest{} in birdie() since they generally take more than 5
  seconds to run. 
