
# check types
check_vec = function(x) (is.character(x) | is.factor(x)) && !any(is.na(x))

# Dev-facing helper for viewing
print_cond = function(x, title=NULL) {
    p_r = parent.frame()$p_r
    out = 100 * x %*% diag(1 / p_r)
    colnames(out) = names(p_r)
    if (!is.null(title)) cat(toupper(title), "\n")
    print(round(out))
}


eval_log_score = function(..., eps=1e-9) {
    log_score = function(m, eps=1e-9) {
        R_vec = as.integer(voters$race)
        mean(log(map_dbl(seq_len(nrow(voters)), \(i) m[i, R_vec[i]] + eps)))
    }

    scores = rlang::list2(...) %>%
        vapply(log_score, numeric(1), eps) %>%
        sort(decreasing=TRUE)
    tibble(method = names(scores),
           log_score = scores)
}
