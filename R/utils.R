check_convergence <- function(est, last, abstol, reltol) {
    diff = abs(est - last)
    absest = abs(est)
    ok = absest > 0

    (max(diff) < abstol) || (max(diff[ok] / absest[ok]) < reltol)
}

# shim. Only works if remaining columns all uniquely identify
pivot_wider_tiny <- function(x, names_from="name", values_from="value") {
    form = as.formula(paste0(values_from, "~", names_from))
    d_id = distinct(select(x, -all_of(c(names_from, values_from))))
    if (nrow(x) == n_distinct(x[[names_from]])) {
        as_tibble(cbind(d_id, unstack(rbind(x, x), form)[1, ]))
    } else {
        as_tibble(cbind(d_id, unstack(x, form)))
    }
}

# converts data frame or vector to integer vector of grouping indices
# e.g. c('a', 'a', 'b', 'a', 'c', 'b') becomes c(1, 1, 2, 1, 3, 2)
to_unique_ids = function(x) {
    as.integer(as.factor(vctrs::vec_duplicate_id(x)))
}


# check types
check_vec = function(x) (is.character(x) | is.factor(x)) && !any(is.na(x))

# Dev-facing helper for viewing
to_cond = function(x) {
    p_r = parent.frame()$p_r
    out = x %*% diag(1 / colSums(x))
    colnames(out) = names(p_r)
    rownames(out) = rownames(x)
    out
}

# Dev-facing helper for viewing
print_cond = function(x, title=NULL, digits=1) {
    p_r = parent.frame()$p_r
    out = 100 * x %*% diag(1 / colSums(x))
    colnames(out) = names(p_r)
    rownames(out) = rownames(x)
    if (!is.null(title)) cat(toupper(title), "\n")
    print(round(out, digits))
}
