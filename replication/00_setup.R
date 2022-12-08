suppressMessages({
    library(here)
    devtools::load_all(here())
    library(tidyverse)
    library(scales)
    library(wacolors)
    library(patchwork)
    library(geomtextpath)
})

races = c(white="White", black="Black", hisp="Hispanic", asian="Asian",
          aian="Native", other="Other")

PAL_R = wacolors$rainier
names(PAL_R) = NULL

theme_paper = function() {
    theme_bw(base_family="Times", base_size=10) +
        theme(plot.background=element_blank())
}



tidy_thresh = function(r_probs, x, nm=deparse(substitute(x)), prefix="pr_") {
    x = as.factor(x)
    N = length(x)

    r_names =  substring(colnames(r_probs), nchar(prefix)+1L)
    stopifnot(nrow(r_probs) == N)

    r_est = predict(r_probs)
    m = prop.table(table(x, r_est), 2)

    out = tibble(X = rep(levels(x), ncol(m)),
                 race = rep(r_names, each=nrow(m)),
                 estimate = as.numeric(m))
    names(out)[1] = nm
    out
}

tidy_ols = function(r_probs, x, nm=deparse(substitute(x)), prefix="pr_") {
    x = as.factor(x)
    N = length(x)

    r_names = substring(colnames(r_probs), nchar(prefix)+1L)
    r_probs = as.matrix(select(r_probs, starts_with(prefix)))
    stopifnot(nrow(r_probs) == N)

    m = lapply(levels(x), function(l) {
        .lm.fit(r_probs, x == l)$coefficients
    })
    m = do.call(rbind, m)

    out = tibble(X = rep(levels(x), ncol(m)),
                 race = rep(r_names, each=nrow(m)),
                 estimate = as.numeric(m))
    names(out)[1] = nm
    out
}

tidy_true = function(race, x, nm=deparse(substitute(x))) {
    x = as.factor(x)
    N = length(x)

    r_names = levels(race)
    stopifnot(length(race) == N)

    m = prop.table(table(x, race), 2)

    out = tibble(X = rep(levels(x), ncol(m)),
                 race = rep(r_names, each=nrow(m)),
                 estimate = as.numeric(m))
    names(out)[1] = nm
    out
}

plot_calib = function(r_probs, race, group="white", bins=16) {
    y_fit = r_probs[[str_c("pr_", group)]]
    tibble(fitted = cut(y_fit, bins),
           y_res = race) %>%
        group_by(fitted) %>%
        summarize(mean = mean(y_res == group),
                  n = n()) %>%
        mutate(fitted = parse_number(as.character(fitted))) %>%
        ggplot(aes(x=fitted, y=mean, size=n)) +
        geom_point() +
        labs(title=str_c("Race: ", group)) +
        scale_size_area(guide="none")  +
        theme_bw()
}

# Helpful function adapted from: https://stat.ethz.ch/pipermail/r-help/2005-September/079872.html
#' Calculate AUC
#'
#' @param x the predictor
#' @param y a binary indicator
#'
#' @returns the scalar AUC
fastAUC <- function(x, y) {
    x1 = x[y==1]; n1 = length(x1);
    x2 = x[y==0]; n2 = length(x2);
    r = rank(c(x1, x2))
    exp(log(sum(r[1:n1]) - n1*(n1 + 1)/2) - log(n1) - log(n2))
}

calc_roc <- function(d_pr) {
    imap(p_r, function(x, r) {
        fastAUC(d_pr[[str_c("pr_", r)]], d$race == r)
    }) |>
        as_tibble()
}
