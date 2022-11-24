test_that("Bayes' rule is calculated correctly for BISG", {
    p_sr = matrix(c(0.75, 0.25, 0.9, 0.1), ncol=2)
    p_gxr = matrix(c(0.5, 0.5, 0.25, 0.75), ncol=2)
    p_r = c(0.5, 0.5)

    num = p_sr[c(1, 1), ] * p_gxr[c(1, 2), ] %*% diag(p_r)
    correct = num / rowSums(num)
    expect_equal(correct, calc_bayes_bisg(c(1, 1), c(1, 2), p_sr, p_gxr, p_r))

    num = p_sr[c(2, 1), ] * p_gxr[c(1, 2), ] %*% diag(p_r)
    correct = num / rowSums(num)
    expect_equal(correct, calc_bayes_bisg(c(2, 1), c(1, 2), p_sr, p_gxr, p_r))
})

test_that("Bayes' rule is calculated correctly in general", {
    p_yr = matrix(c(0.75, 0.25, 0.1, 0.9), ncol=2)
    prior = matrix(c(seq(0, 1, 0.1), seq(1, 0, -0.1)), ncol=2)

    y = c(rep(1:2, 5), 1)
    num = p_yr[y, ] * prior
    correct = num / rowSums(num)
    expect_equal(correct, calc_bayes(y, rep_along(y, 1), as.numeric(t(p_yr)), prior, 1, 2))
})
