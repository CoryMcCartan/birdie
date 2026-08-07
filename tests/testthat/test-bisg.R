log_score <- function(x, R) {
    pr_act = as.matrix(x)[cbind(1:nrow(x), as.integer(R))]
    pr_act[pr_act == 0] = 1e-6
    mean(log(pr_act))
}

test_that("BISG checks input", {
    data(pseudo_vf)

    expect_s3_class(bisg(~ nm(last_name), pseudo_vf), "bisg")
    expect_s3_class(bisg(~ nm(last_name) + zip(zip), pseudo_vf), "bisg")
    expect_error(expect_warning(
        bisg(~ nm(last_name) + zip(zip) + race, pseudo_vf)
        ))
    expect_error(bisg(~ nm(last_name) + zip(zip) + state(zip), pseudo_vf))
    expect_error(bisg(~ nm(last_name) + zip(zip) + state(race), pseudo_vf))
    expect_error(bisg(~ nm(last_name) + zip + zip, pseudo_vf))
    expect_error(bisg(~ zip(zip), pseudo_vf))
})

test_that("BISG is robust to ordering", {
    data(pseudo_vf)
    p_r = p_r_natl(2010)

    expect_s3_class(bisg(~ nm(last_name), pseudo_vf, p_r=rev(p_r)), "bisg")
    expect_error(bisg(~ nm(last_name), pseudo_vf, p_r=p_r[1:3]), "doesn't match")
    expect_error(bisg(~ nm(last_name), pseudo_vf, p_r=unname(p_r)), "must have names")
})

test_that("National race demographics cover the documented years", {
    expect_equal(p_r_natl(), p_r_natl(2024))
    expect_equal(sum(p_r_natl(2005)), 1, tolerance=1e-6)
    expect_equal(sum(p_r_natl(2024)), 1, tolerance=1e-6)
    expect_error(p_r_natl(2004), "between 2005 and 2024")
    expect_error(p_r_natl(2025), "between 2005 and 2024")
})


test_that("BISG results are directionally correct", {
    d = tibble(S=rep(c("Hernandez", "Mc Cartan"), 2),
               G=rep(c("78501", "50112"), each=2))
    p_r = p_r_natl(2010)
    res = bisg(~ nm(S) + zip(G), data=d, p_r=p_r)

    expect_gt(res$pr_hisp[1], 0.98)
    expect_gt(res$pr_white[4], 0.98)
    expect_gt(res$pr_hisp[1], res$pr_hisp[3])
    expect_gt(res$pr_white[4], res$pr_white[2])
    expect_gt(res$pr_hisp[1], res$pr_hisp[2])
    expect_gt(res$pr_white[4], res$pr_white[3])

    expect_gt(bisg(~ nm(S) + zip(G), data=tibble(S="LOCKLEAR", G="28715"), p_r=p_r)$pr_aian, 0.4)
    expect_gt(bisg(~ nm(S) + zip(G), data=tibble(S="WASHINGTON", G="98118"), p_r=p_r)$pr_black, 0.9)
    expect_gt(bisg(~ nm(S) + zip(G), data=tibble(S="XU", G="98118"), p_r=p_r)$pr_asian, 0.95)
})

test_that("BISG results do not depend on the prediction batch", {
    p_r = p_r_natl(2010)
    d_one = tibble(S="HERNANDEZ", G="78501")
    d_many = tibble(S=c("HERNANDEZ", rep("SMITH", 20)),
                    G=c("78501", rep("50112", 20)))

    b_one = bisg(~ nm(S) + zip(G), data=d_one, p_r=p_r)
    b_many = bisg(~ nm(S) + zip(G), data=d_many, p_r=p_r)

    expect_equal(as.numeric(b_one[1, ]), as.numeric(b_many[1, ]))
})

test_that("BISG uses fixed reference likelihoods with custom tables", {
    d = tibble(S=c("A", "B"), G=c("x", "y"))
    p_rs = tibble(S=c("A", "B"), r1=c(60, 40), r2=c(20, 80))
    p_rgx = tibble(G=c("x", "y"), r1=c(30, 70), r2=c(80, 20))
    p_r = c(r1=0.4, r2=0.6)

    out = bisg(~ nm(S) + G, d, p_r=p_r, p_rs=p_rs, p_rgx=p_rgx)
    p_sr = as.matrix(p_rs[, -1]) / colSums(as.matrix(p_rs[, -1]))
    p_rg = as.matrix(p_rgx[, -1]) / rowSums(as.matrix(p_rgx[, -1]))
    q_r = colSums(as.matrix(p_rgx[, -1])) / sum(as.matrix(p_rgx[, -1]))
    expected = p_sr * sweep(p_rg, 2, q_r, "/") * rep(p_r, each=2)
    expected = expected / rowSums(expected)
    colnames(expected) = paste0("pr_", colnames(expected))

    expect_equal(as.matrix(out), expected)
    expect_equal(attr(out, "q_r"), q_r)

    vars = parse_bisg_form(~ nm(S) + G, d)
    expect_equal(make_name_tbl_vec(vars, p_r, p_rs, for_me=TRUE)$p_sr,
                 as.matrix(p_rs[, -1]))
    expect_equal(unname(make_gx_tbl_vec(vars, p_r, p_rgx)$p_rgx),
                 unname(as.matrix(p_rgx[, -1])))
})

test_that("BISG estimates the target race prior by EM", {
    d = tibble(S=rep(c("A", "B"), c(3, 7)), G=rep(c("x", "y"), c(3, 7)))
    p_rs = tibble(S=c("A", "B"), r1=c(90, 10), r2=c(10, 90))
    p_rgx = tibble(G=c("x", "y"), r1=c(80, 20), r2=c(20, 80))

    out = bisg(~ nm(S) + G, d, p_r="estimate", p_rs=p_rs, p_rgx=p_rgx)

    expect_equal(unname(attr(out, "p_r")), unname(colMeans(out)), tolerance=1e-8)
})

test_that("Unmatched geography does not update the target prior", {
    d = tibble(S="A", G="missing")
    p_rs = tibble(S="A", r1=1, r2=1)
    p_rgx = tibble(G=c("x", "y"), r1=c(40, 10), r2=c(10, 40))
    p_r = c(r1=0.8, r2=0.2)

    out = bisg(~ nm(S) + G, d, p_r=p_r, p_rs=p_rs, p_rgx=p_rgx)

    expect_equal(as.numeric(out[1, ]), unname(p_r))
})

test_that("BISG works with custom p_rs", {
    data(pseudo_vf)

    p_rs = census_surname_table(proc_name(pseudo_vf$last_name), "last_name", count=TRUE)

    b0 = bisg(~ nm(last_name) + zip(zip), data=pseudo_vf)
    b1 = bisg(~ nm(last_name) + zip(zip), data=pseudo_vf, p_rs=p_rs)

    expect_lt(mean(abs(colMeans(b1) - colMeans(b0))), 0.05)
    expect_s3_class(b0, "bisg")
    expect_s3_class(b1, "bisg")
})

test_that("Measurement error BISG model works", {
    data("pseudo_vf")

    p_r = prop.table(table(pseudo_vf$race))

    pr_0 = bisg(~ nm(last_name) + zip(zip), data=pseudo_vf, p_r=p_r)
    pr_me = bisg_me(~ nm(last_name) + zip(zip), data=pseudo_vf, p_r=p_r,
                    warmup=100, iter=1000)

    expect_true(all(diag(cor(pr_0, pr_me))[1:5] > 0.9))

    # ME within 2% or better
    expect_gt(log_score(pr_me, pseudo_vf$race), log_score(pr_0, pseudo_vf$race) - 0.02)

    # should work with estimate
    expect_no_error(
        bisg_me(~ nm(last_name) + zip(zip), data=pseudo_vf, p_r="estimate",
                warmup=10, iter=100)
    )
})

test_that("BISG results match `wru`", {
    skip_on_cran()
    skip_on_ci()
    skip_if_offline("github.com")
    skip_if_not_installed("wru")

    data("pseudo_vf")
    pseudo_vf$dummy = factor(1)

    p_r = c(white=0.630, black=0.121, hisp=0.173,
          asian=0.0478, aian=0.0072, other=0.0210)
    p_r = p_r / sum(p_r)

    m_wru = pseudo_vf |>
        dplyr::rename(surname=last_name) |>
        dplyr::mutate(state="NC") |>
        wru::predict_race(surname.only=TRUE) |>
        dplyr::select("pred.whi":"pred.oth") |>
        suppressMessages() |>
        as.matrix()

    m_birdie = bisg(~ nm(last_name), data=pseudo_vf, p_r=p_r) |>
        dplyr::select(-pr_aian) |>
        as.matrix()
    m_birdie = m_birdie / rowSums(m_birdie)

    expect_true(all(diag(cor(m_birdie, m_wru)[1:4, 1:4]) >= 0.9))
})
