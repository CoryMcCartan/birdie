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
#' @param object,x A `birdie` model object
#' @param data A data frame to augment with `Pr(R | Y, X, S)` probabilities
#' @param ... Potentially further arguments passed from other methods
#'
#' @return Varies, depending on the method. See generic functions' documentation
#'   for details.
#'
#' @examples
#' methods(class="birdie")
#'
#' @concept estimators
#' @name birdie-class
NULL


#' @param subgroup If `TRUE`, return subgroup-level (rather than marginal)
#'   coefficient estimates as a 3D array.
#'
#' @describeIn birdie-class Return estimated outcome-given-race distributions.
#' @export
coef.birdie <- function(object, subgroup=FALSE, ...) {
    if (isFALSE(subgroup)) {
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

#' @describeIn birdie-class Return the residuals for the outcome variable.
#'   Useful in sensitivity analyses and to get an idea of how well race,
#'   location, names, etc. predict the outcome.
#' @param x_only if `TRUE`, calculate fitted values using covariates only (i.e.,
#'   without using surnames).
#' @export
residuals.birdie <- function(object, x_only=FALSE, ...) {
    m = as.numeric(coef.birdie(object, subgroup=TRUE))
    r_probs = if (isFALSE(x_only)) {
        as.matrix(object$p_ryxs)
    } else {
        tmp = attr(object$p_ryxs, "p_rgx")
        if (is.null(tmp)) {
            cli_abort(c("Cannot generate residuals with {.arg x_only=TRUE}.",
                        "i"="Missing Pr(R | G, X) table.",
                        ">"="Generate your {.fn bisg} predictions with
                        {.arg save_rgx=TRUE}."))
        }
        as.matrix(tmp[attr(object$p_ryxs, "gx"), ])
    }
    n_y = nlevels(object$y)
    y = as.integer(object$y)

    out = do.call(cbind, lapply(seq_len(n_y), function(i) {
        (y == i) - resid_mult(m, object$vec_gx, r_probs, i, n_y)
    }))
    colnames(out) = levels(object$y)
    rownames(out) = NULL
    out
}

#' @describeIn birdie-class Create point predictions of individual race. Returns
#'   factor vector of individual race labels. Strongly not recommended for any
#'   kind of inferential purpose, as biases may be extreme and in unpredictable
#'   directions.
#' @inheritParams predict.bisg
#' @export
predict.birdie <- function(object, adj=NULL, ...) {
    predict.bisg(object$p_ryxs, adj=adj, ...)
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


#' @describeIn birdie-class Visualize the estimated conditional distributions
#'   for a BIRDiE model.
#' @param log If `TRUE`, plot estimated probabilities on a log scale.
#' @method plot birdie
#' @export
plot.birdie <- function(x, log=FALSE, ...) {
    m = coef.birdie(x)

    resp = x$y_name
    ylab = str_c("Pr(", resp, " | Race)")
    main = paste("Estimates of", resp, "by race")
    n_y = nrow(m)

    PAL_RAINIER = c("#465177", "#E4C22B", "#965127", "#29483A", "#759C44",
                             "#9FB6DA", "#DF3383")
    PAL_PUGET = c("#1D3024", "#123B2D", "#00473E", "#005252", "#005B66",
                "#00657D", "#386B91", "#5F6EA3", "#8172B1", "#A074B8",
                "#B87DB8", "#CB85B6", "#D992B2", "#E59FAE", "#EEADAB")

    if (n_y <= 7) {
        pal = PAL_RAINIER[seq_len(n_y)]
    } else {
        pal = colorRampPalette(PAL_PUGET)(nrow(m))
    }

    barplot(m, names.arg=toupper(colnames(m)), log=if (log) "y" else "",
         cex.names=0.85, border=NA, space=c(0.1, 0.9), axis.lty=0,
         ylab=ylab, xlab="Race", col=pal, beside=TRUE,
         args.legend=list(x="topright", box.lwd=0, cex=0.9, y.intersp=0.1,
                          horiz=TRUE, inset=c(0, -0.05)),
         legend.text=TRUE, ...)

    invisible(m)
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
tidy.birdie <- function(x, subgroup=FALSE, ...) {
    d = dim(x$map_sub)
    if (isFALSE(subgroup) || d[3] == 1) {
        m = x$map
        out = tibble(X = rep(rownames(m), ncol(m)),
                     race = rep(colnames(m), each=nrow(m)),
                     estimate = as.numeric(m))
        names(out)[1] = x$y_name
    } else {
        m = x$map_sub
        dn = dimnames(m)
        out = tibble(X = rep(dn[[1]], d[2]*d[3]),
                     race = rep(dn[[2]], d[3], each=d[1]),
                     .i = rep(seq_len(d[3]), each=d[1]*d[2]),
                     estimate = as.numeric(m))
        names(out)[1] = x$y_name

        gx = x$tbl_gx
        gx$.i = seq_len(d[3])
        out = left_join(gx, out, by=".i")
        out$.i = NULL
    }

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
           entropy.pre.med = median(x$entropy$pre),
           entropy.post.med = median(x$entropy$post),
           nobs = nrow(m),
           ngrp = dim(x$map_sub)[3])
}

#' @describeIn birdie-class Augment data with individual race predictions from a BIRDiE model.
#' @method augment birdie
#' @export
augment.birdie <- function(x, data, ...) {
    bind_cols(data, fitted.birdie(x))
}

#' @describeIn birdie-class Extract the formula used to specify a BIRDiE model.
#' @method formula birdie
#' @export
formula.birdie <- function(x, ...) {
    NextMethod()
}

#' @describeIn birdie-class Return the BIRDiE complete-data model family.
#' @method nobs birdie
#' @export
family.birdie <- function(object, ...) {
    object$family
}

#' @describeIn birdie-class Return the number of observations used to fit a BIRDiE model.
#' @method nobs birdie
#' @export
nobs.birdie <- function(object, ...) {
    nrow(object$p_ryxs)
}

#' @describeIn birdie-class Return the estimated variance-covariance matrix for
#'   the BIRDiE model estimates, if available.
#' @method vcov birdie
#' @export
vcov.birdie <- function(object, ...) {
    object$vcov
}



comma <- function(x) format(x, big.mark=',')

#' @describeIn birdie-class Print a summary of the model fit.
#' @export
print.birdie <- function(x, ...) {
    cli::cli_text(switch(
        x$algo$model,
        cat_dir = "Categorical-Dirichlet {.pkg BIRDiE} model",
        cat_mixed = "Categorical mixed-effects {.pkg BIRDiE} model",
        "{.pkg BIRDiE} model"
    ))
    cli::cat_line("Formula: ", deparse(x$call$formula))
    cli::cat_line("   Data: ", deparse(x$call$data))
    cli::cli_text("Number of obs: {comma(x$N)};
                  groups: {comma(dim(x$map_sub)[3])}")

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

#' @describeIn birdie-class Print a more detailed summary of the model fit.
#' @export
summary.birdie <- function(object, ...) {
    cli::cli_text(switch(
        object$algo$model,
        cat_dir = "Categorical-Dirichlet {.pkg BIRDiE} model",
        cat_mixed = "Categorical mixed-effects {.pkg BIRDiE} model",
        "{.pkg BIRDiE} model"
    ))
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
    cli::cli_text("Number of groups: {comma(dim(object$map_sub)[3])}")
    cat("\n")

    p_r = colMeans(object$p_ryxs)
    cli::cat_line("Entropy decrease from marginal race distribution:")
    print(summary(entropy(p_r) - object$entropy$post))
    cli::cat_line("Entropy decrease from BISG probabilities:")
    print(summary(object$entropy$pre - object$entropy$post))
    cat("\n")

    cli::cli_text("Estimated outcome-by-race distribution:")
    m = round(object$map, 3)
    print(m)

    invisible(object)
}

# WEIGHTED ESTIMATOR PRINT/SUMMARY GENERICS ARE IN `weighted.R`


#' @describeIn bisg Summarize predicted race probabilities. Returns vector of individual entropies.
#' @param object An object of class `bisg`, the result of running [bisg()].
#' @param ... Additional arguments to generic methods (ignored).
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

#' @describeIn bisg Create point predictions of individual race. Returns factor
#'   vector of individual race labels. Strongly not recommended for any kind of
#'   inferential purpose, as biases may be extreme and in unpredictable
#'   directions.
#' @param adj A point in the simplex that describes how BISG probabilities
#'   will be thresholded to produce point predictions. The probabilities are
#'   divided by `adj`, then the racial category with the highest probability is
#'   predicted. Can be used to trade off types of prediction error. Must be
#'   nonnegative but will be normalized to sum to 1. The default is to make no
#'   adjustment.
#' @export
predict.bisg <- function(object, adj=NULL, ...) {
    n_r = ncol(object)
    races = stringr::str_sub(colnames(object), 4)

    if (is.null(adj)) {
        adj = rep(1, n_r)
    } else {
        if (any(adj < 0)) {
            cli::cli_abort("{.arg thresh} must be nonnegative.")
        }
    }
    adj = adj / sum(adj)

    m = as.matrix(object) / rep(adj, each=nrow(object))

    factor(max.col(m), levels=seq_len(n_r), labels=races)
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
    attr(new, "p_rgx") = attr(old, "p_rgx")
    attr(new, "gx") = attr(old, "gx")
    attr(new, "method") = "birdie"
    class(new) = c("bisg", class(new))
    new
}
