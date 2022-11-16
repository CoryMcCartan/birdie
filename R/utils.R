
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
