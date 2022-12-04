#' Class "birdie" of BIRDiE Models
#'
#' The output of [birdie()] is an object of class `birdie`, which supports
#' many generic functions. Notably `coef.birdie()` returns the main model
#' estimates of outcome given race, and `fitted.birdie()` returns a table
#' analogous to the output of [bisg()] with updated race probabilities.
#'
#' The internal structure of `birdie` objects is not designed to be accessed
#' directly. The generics listed here should be used instead.
#'
#' @param object,x A `birdie` object
#' @param data A data frame to agument with `Pr(R | Y, X, S)` probabilities
#' @param ... Potentially further arguments passed from other methods
#'
#' @name birdie-class
NULL


#' @param complete If `TRUE`, return group-level (rather than marginal)
#'   coefficient estimates as a 3D array.
#'
#' @describeIn birdie-class Return estimated outcome-given-race distributions.
#' @export
coef.birdie <- function(object, complete=FALSE, ...) {
    if (isFALSE(complete)) {
        object$map
    } else {
        object$map_sub
    }
}

#' @describeIn birdie-class Return an updated race probability table. [bisg()]
#'   estimates `Pr(R | G, X, S)`; this table estimates `Pr(R | Y, G, X, S)`.
#' @export
fitted.birdie <- function(object, ...) {
    object$p_ryxs
}

#' @describeIn birdie-class Alias for `fitted.birdie`.
#' @export
predict.birdie <- function(object, ...) {
    fitted.birdie(object)
}

#' @param nsim The number of vectors to simulate. Defaults to 1.
#' @param seed Used to seed the random number generator. See [stats::simulate()].
#'
#' @describeIn birdie-class Simulate race from the posterior distribution
#'   `Pr(R | Y, G, X, S, Theta-MAP)`. Does not account for uncertainty in model
#'   parameters.
#' @export
simulate.birdie <- function(object, nsim = 1, seed = NULL, ...) {
    if (!is.null(seed)) set.seed(seed)

    m = as.matrix(object$p_ryxs)

    out = do.call(cbind, lapply(seq_len(nsim), function(i) {
        mat_rcatp(m)
    }))

    r_names = substring(colnames(m), nchar(object$prefix)+1L)
    out = factor(out, levels=seq_len(ncol(m)), labels=r_names)
    dim(out) = c(nrow(m), nsim)

    out
}

#' @importFrom generics tidy glance augment
#' @export
generics::tidy
#' @export
generics::glance
#' @export
generics::augment

#' @describeIn birdie-class Put BIRDiE model coefficients in a tidy format.
#' @method tidy birdie
#' @export
tidy.birdie <- function(x, ...) {
    m = x$map

    out = tibble(X = rep(rownames(m), ncol(m)),
                 race = rep(colnames(m), each=nrow(m)),
                 estimate = as.numeric(m))
    names(out)[1] = rlang::expr_name(rlang::f_lhs(x$call$formula))

    out
}

#' @describeIn birdie-class Glance at a BIRDiE model.
#' @method glance birdie
#' @export
glance.birdie <- function(x, ...) {
    m = x$p_ryxs

    p_r = colMeans(m)
    ents = entropy(m)

    tibble(entropy.marg = entropy(p_r),
           entropy.post.med = median(entropy(m)),
           nobs = nrow(m),
           ngrp = dim(x$map_sub)[1])
}

#' @describeIn birdie-class Augment data with individual race predictions from a BIRDiE model.
#' @method augment birdie
#' @export
augment.birdie <- function(x, data, ...) {
    bind_cols(data, fitted.birdie(x))
}


comma <- function(x) format(x, big.mark=',')

#' @describeIn birdie-class Print a summary of the model fit.
#' @export
print.birdie <- function(x, ...) {
    cli::cli_text("{.pkg BIRDiE} model fit with method = {.arg {x$algo$method}}")
    cli::cat_line("Formula: ", deparse(x$call$formula))
    cli::cat_line("   Data: ", deparse(x$call$data))
    cli::cli_text("Number of obs: {comma(x$N)};
                  groups: {comma(dim(x$map_sub)[1])}")

    cli::cli_text("Estimated distribution:")
    m = round(x$map, 3)
    print(m)
}

entropy <- function(x) {
    if (is.data.frame(x)) {
        x = as.matrix(x)
    } else if (!is.matrix(x)) {
        x = matrix(x, nrow=1)
    }
    -rowSums(x * log(x), na.rm=TRUE)
}

#' @describeIn birdie-class Print a summary of the model fit.
#' @export
summary.birdie <- function(object, ...) {
    cli::cli_text("{.pkg BIRDiE} model fit with method = {.arg {object$algo$method}}")
    cli::cat_line("Formula: ", deparse(object$call$formula))
    cli::cat_line("   Data: ", deparse(object$call$data))
    cat("\n")

    if (object$algo$converge) {
        t0 = Sys.time()
        tm = format(difftime(t0 + object$algo$runtime, t0), digits=2)
        cli::cli_text("{object$algo$iters} iterations and {tm} to convergence")
    } else {
        cli::cli_alert_danger("Model did not converge")
    }
    cat("\n")

    cli::cli_text("Number of observations: {comma(object$N)}")
    cli::cli_text("Number of groups: {comma(dim(object$map_sub)[1])}")
    cat("\n")

    p_r = colMeans(object$p_ryxs)
    cli::cat_line("Entropy decrease from marginal race distribution:")
    print(summary(entropy(p_r) - entropy(object$p_ryxs)))
    cat("\n")

    cli::cli_text("Estimated outcome-by-race distribution:")
    m = round(object$map, 3)
    print(m)

    invisible(object)
}

#' @describeIn bisg Summarize predicted race probabilities. Returns vector of individual entropies.
#' @param object An object of class `bisg`, the result of running [bisg()].
#' @param ... Additional arguments to `summary()` (ignored).
#' @export
summary.bisg <- function(object, p_r=NULL, ...) {
    cli::cli_text("{.pkg BISG} individual race probabilities")
    if (inherits(object, "bisg_me")) {
        cli::cat_line("Estimated with measurement error model")
    }
    cat("\n")

    if (is.null(p_r)) {
        cli::cat_line("Implied marginal race distribution:")
        p_r = colMeans(object)
    } else {
        cli::cat_line("Given marginal race distribution:")
    }
    print(round(p_r, 3))
    cat("\n")

    cli::cat_line("Entropy decrease from marginal distribution:")
    ents <- entropy(object)
    print(summary(entropy(p_r) - ents))

    invisible(ents)
}

reconstruct <- function(new, old, ...) {
    UseMethod("reconstruct")
}

reconstruct.default <- function(new, old, ...) {
    new
}

reconstruct.bisg <- function(new, old, cat_nms=NULL, ...) {
    attr(new, "S_name") = attr(old, "S_name")
    attr(new, "GX_names") = c(cat_nms, attr(old, "GX_names"))
    attr(new, "p_r") = attr(old, "p_r")
    attr(new, "method") = "birdie"
    class(new) = c("bisg", class(new))
    new
}
