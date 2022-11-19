# shim. Only works if remaining columns all uniquely identify
pivot_wider_tiny <- function(x, names_from="name", values_from="value") {
    form = as.formula(paste0(values_from, "~", names_from))
    d_id = distinct(select(x, c(-names_from, -values_from)))
    as_tibble(cbind(d_id, unstack(x, form)))
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
