data(pseudo_vf)
r_probs = bisg(~ nm(last_name) + zip(zip), data=pseudo_vf)
zip_patched = proc_zip(pseudo_vf$zip)

test_that("weighted estimator runs", {
    expect_s3_class(
        est_weighted(r_probs, turnout ~ 1, data=pseudo_vf),
        "est_weighted")

    expect_s3_class(
        est_weighted(r_probs, turnout ~ proc_zip(zip), data=pseudo_vf),
        "est_weighted")
})

test_that("weighted estimator is calculated correctly", {
    est = coef(est_weighted(r_probs, turnout ~ 1, data=pseudo_vf))

    x = pseudo_vf$turnout
    est0 = do.call(rbind, lapply(levels(x), function(l) {
        colSums(r_probs * (x == l)) / colSums(r_probs)
    }))
    dimnames(est0) = dimnames(est)

    expect_equal(est, est0)
})
