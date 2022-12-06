data(pseudo_vf)
r_probs = bisg(~ nm(last_name) + zip(zip), data=pseudo_vf)
zip_patched = proc_zip(pseudo_vf$zip)
X = sample(c("a", "b", "c"), nrow(pseudo_vf), replace=TRUE)

test_that("BIRDiE models fit", {
    ctrl = birdie.ctrl(abstol=1e-3, max_iter=200)

    expect_s3_class(
        birdie(r_probs, turnout ~ 1, data=pseudo_vf, ctrl=ctrl),
        "birdie")

    expect_s3_class(
        birdie(r_probs, turnout ~ proc_zip(zip), data=pseudo_vf, ctrl=ctrl),
        "birdie")

    expect_s3_class(
        birdie(r_probs, turnout ~ (1 | proc_zip(zip)), data=pseudo_vf, ctrl=ctrl),
        "birdie")

    expect_s3_class(
        suppressWarnings(
        birdie(r_probs, turnout ~ X, model="mmm", data=pseudo_vf, ctrl=ctrl),
        ), "birdie")
})

test_that("BIRDiE produces correct and stable output", {
    skip_on_cran()

    out = birdie(r_probs, turnout ~ 1, data=pseudo_vf)
    p_r =  with(pseudo_vf, prop.table(table(race)))
    p_xr = with(pseudo_vf, prop.table(table(turnout, race), 2))

    tv = sum(abs(coef(out) - p_xr) %*% diag(p_r)) / 2

    expect_lt(tv, 0.015)

    expect_snapshot_value(coef(out), style="json2", tolerance=1e-5)
})

test_that("BIRDiE catches errors", {
    ctrl = birdie.ctrl(abstol=1e-3, max_iter=200)

    expect_error(
        birdie(r_probs, turnout ~ 1, data=pseudo_vf, prior=list(), ctrl=ctrl),
        "alpha"
    )

    expect_error(
        birdie(r_probs, turnout ~ 1, data=pseudo_vf,
               prior=list(alpha=matrix(1, nrow=3, ncol=6)), ctrl=ctrl),
        "rows"
    )

    expect_warning(
        birdie(r_probs, turnout ~ zip_patched, data=pseudo_vf, ctrl=ctrl),
        "Missing"
    )

    attr(r_probs, "GX_names") = character(0) # hide the "missing" warning for next test
    expect_warning(
        birdie(r_probs, turnout ~ last_name, data=pseudo_vf, ctrl=ctrl),
        "Last name"
    )
})
