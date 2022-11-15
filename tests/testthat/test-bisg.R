log_score <- function(x, R) {
    pr_act = as.matrix(x)[cbind(1:nrow(x), as.integer(R))]
    pr_act[pr_act == 0] = 1e-6
    mean(log(pr_act))
}


test_that("BISG results are directionally correct", {
    d = tibble(S=rep(c("Hernandez", "Mc Cartan"), 2),
               G=rep(c("78501", "50112"), each=2))
    res = predict_race_sgz(S, G, data=d)

    expect_gt(res$pr_hisp[1], 0.99)
    expect_gt(res$pr_white[4], 0.99)
    expect_gt(res$pr_hisp[1], res$pr_hisp[3])
    expect_gt(res$pr_white[4], res$pr_white[2])
    expect_gt(res$pr_hisp[1], res$pr_hisp[2])
    expect_gt(res$pr_white[4], res$pr_white[3])

    expect_gt(predict_race_sgz(S, G, data=tibble(S="LOCKLEAR", G="28715"))$pr_aian, 0.5)
    expect_gt(predict_race_sgz(S, G, data=tibble(S="WASHINGTON", G="98118"))$pr_black, 0.9)
    expect_gt(predict_race_sgz(S, G, data=tibble(S="CHEN", G="98118"))$pr_asian, 0.95)
})

test_that("Measurement error BISG model works", {
    data("pseudo_vf")

    p_r = prop.table(table(pseudo_vf$race))

    pr_0 = predict_race_sgz(last_name, zip, data=pseudo_vf, p_r=p_r)
    pr_me = predict_race_sgz_me(last_name, zip, data=pseudo_vf, p_r=p_r,
                                warmup=100, iter=1000)

    expect_true(all(diag(cor(pr_0, pr_me))[1:5] > 0.9))

    # ME within 2% or better
    expect_gt(log_score(pr_me, pseudo_vf$race), log_score(pr_0, pseudo_vf$race) - 0.02)
})

test_that("BISG results match `wru`", {
    data("pseudo_vf")
    pseudo_vf$dummy = factor(1)

    p_r = c(white=0.630, black=0.121, hisp=0.173,
          asian=0.0478, aian=0.0072, other=0.0210)
    p_r = p_r / sum(p_r)

    m_wru = pseudo_vf |>
        dplyr::rename(surname=last_name) |>
        dplyr::mutate(state="NC") |>
        wru::predict_race(surname.only=TRUE) |>
        dplyr::select(.data$pred.whi:.data$pred.oth) |>
        as.matrix()

    p_rgz = matrix(p_r, nrow=1) |>
        `colnames<-`(names(p_r)) |>
        as_tibble() |>
        dplyr::mutate(dummy = "1", .before="white")

    m_birdie = predict_race_sgz(last_name, G=dummy, data=pseudo_vf, p_rgz=p_rgz) |>
        dplyr::select(-pr_aian) |>
        as.matrix()
    m_birdie = m_birdie / rowSums(m_birdie)

    expect_true(all(diag(cor(m_birdie, m_wru)[1:4, 1:4]) >= 0.9))
})
