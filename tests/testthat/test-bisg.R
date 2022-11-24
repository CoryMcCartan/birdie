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

    expect_gt(bisg(~ nm(S) + zip(G), data=tibble(S="LOCKLEAR", G="28715"), p_r=p_r)$pr_aian, 0.5)
    expect_gt(bisg(~ nm(S) + zip(G), data=tibble(S="WASHINGTON", G="98118"), p_r=p_r)$pr_black, 0.9)
    expect_gt(bisg(~ nm(S) + zip(G), data=tibble(S="XU", G="98118"), p_r=p_r)$pr_asian, 0.95)
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
})

test_that("BISG results match `wru`", {
    skip_if_offline("github.com")

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
