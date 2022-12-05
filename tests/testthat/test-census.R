skip_on_cran()
skip_if(Sys.getenv("CENSUS_API_KEY") == "")

test_that("Census data can be downloaded", {
    expect_s3_class(census_race_geo_table("us"), "data.frame")
    expect_s3_class(census_race_geo_table("region"), "data.frame")
    expect_s3_class(census_race_geo_table("zcta"), "data.frame")
    expect_s3_class(census_race_geo_table("state"), "data.frame")
    expect_s3_class(census_race_geo_table("puma", year=2020, survey="acs5"),
                    "data.frame")
})
