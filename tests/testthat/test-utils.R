test_that("convergence is checked correctly", {
    expect_false(check_convergence(1, 1.01, 0, 0))
    expect_true(check_convergence(1, 1.01, 0.01001, 0))
    expect_true(check_convergence(1, 1.01, 0, 0.01001))
    expect_false(check_convergence(10, 10.1, 0.01001, 0))
    expect_true(check_convergence(10, 10.1, 0, 0.01001))
})

# NO LONGER USED
# test_that("simplex transformations are correct", {
#     expect_equal(to_simplex(0), rep(0.5, 2))
#     expect_equal(from_simplex(rep(0.5, 2)), 0)
#     expect_equal(to_simplex(rep(0, 3)), rep(0.25, 4))
#     expect_equal(from_simplex(rep(0.25, 4)), rep(0, 3))
#
#     x = rnorm(100)
#     y = to_simplex(x)
#     expect_equal(from_simplex(y), x)
#     expect_equal(to_simplex(from_simplex(y)), y)
# })
