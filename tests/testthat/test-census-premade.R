test_that("2020 Census geography tables are available", {
    race_cols = c("white", "black", "hisp", "asian", "aian", "other")

    state_counts = census_race_geo_table("state", year=2020, counts=TRUE)
    state_props = census_race_geo_table("state", year=2020, counts=FALSE)
    zcta_counts = census_race_geo_table("zcta", year=2020, counts=TRUE)
    zcta_props = census_race_geo_table("zcta", year=2020, counts=FALSE)

    expect_named(state_counts, c("GEOID", race_cols))
    expect_named(zcta_counts, c("GEOID", race_cols))
    expect_equal(nrow(state_counts), 52)
    expect_gt(nrow(zcta_counts), 33000)
    expect_equal(rowSums(state_props[race_cols]), rep(1, nrow(state_props)))
    zcta_totals = rowSums(zcta_counts[race_cols])
    expect_equal(rowSums(zcta_props[zcta_totals > 0, race_cols]),
                 rep(1, sum(zcta_totals > 0)))
    expect_true(all(is.nan(as.matrix(zcta_props[zcta_totals == 0, race_cols]))))
})
