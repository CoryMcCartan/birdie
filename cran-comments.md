## Test environments
* local R installation (macOS 14.5, arm64), R 4.4.0
* macOS 12.6.1 (on GitHub Actions), (release)
* Windows Server 2022 (on GitHub Actions), (release)
* ubuntu 22.04 (on GitHub Actions), (devel)
* ubuntu 22.04 (on GitHub Actions), (release)
* ubuntu 22.04 (on GitHub Actions), (old release R 4.1.3)

## R CMD check results

0 errors | 0 warnings | 2 notes

* This is a new release.

* C++17 is a requirement in order to compile Stan models.

* GNU make is a SystemRequirement in order to compile Stan models as well.

* Examples are \dontrun{} in census_race_geo_table() since they require an API
key and may take some time to run.
