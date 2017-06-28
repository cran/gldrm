#' Confidence intervals for gldrm coefficients
#'
#' Calculates a Wald or likelihood ratio confidence interval for a single gldrm
#' coefficient. Also calculates upper or lower confidence bounds. Wald confidence
#' intervals and bounds are calculated from the standard errors which are available
#' from the gldrm model fit. For likelihood ratio intervals and bounds, a bisection
#' search method is used, which takes longer to run.
#'
#' @param gldrmFit A fitted gldrm model.
#' @param term Character string containing the name of the coefficient of interest.
#' The coefficient names are the names of the beta component of the fitted model
#' object. They can also be obtained from the printed model output. Usually the
#' names match the formula syntax, but can be more complicated for categorical
#' variables and interaction terms.
#' @param test Character string for the type confidence interval. Options are
#' "Wald" or "LRT" (for likelihood ratio).
#' @param level Confidence level of the interval. Should be between zero and one.
#' @param type Character string containing "2-sided" for a two-sided confidence interval,
#' "lb" for a lower bound, or "ub" for an upper bound.
#' @param eps Convergence threshold. Only applies for \code{test = "LRT"}.
#' Convergence is reached when likelihood ratio p-value is within \code{eps} of
#' the target p-value, based on the level of the test. For example, a two-sided
#' 95\% confidence interval has target p-value of 0.025 for both the upper and
#' lower bounds. A 95\% confidence bound has target p-value 0.05.
#' @param maxiter The maximum number of bisection method iterations for likelihood
#' ratio intervals or bounds. For two-sided intervals, \code{maxiter} iterations
#' are allowed for each bound.
#'
#' @return An S3 object of class 'gldrmCI', which is a list of the following items.
#'
#' \itemize{
#' \item \code{term} Coefficient name.
#' \item \code{test} Type of interval or bound - Wald or likelihood ratio.
#' \item \code{level} Confidence level.
#' \item \code{type} Type of interval or bound - two-sided, upper bound, or lower
#' bound.
#' \item \code{cilo}/\code{cihi} Upper and lower interval bounds. One one of the
#' two applies for confidence bounds.
#' \item \code{iterlo}/\code{iterhi} Number of bisection iterations used. Only
#' applies for likelihood ratio intervals and bounds.
#' \item \code{pvallo}/\code{pvalhi} For likelihood ratio intervals and bounds,
#' the p-value at convergence is reported.
#' }
#'
#' @examples
#' data(iris, package="datasets")
#'
#' ### Fit gldrm with all variables
#' lf <- make.link("log")
#' linkfun <- lf$linkfun  # this is equivalent to function(mu) log(mu)
#' linkinv <- lf$linkinv  # this is equivalent to function(eta) exp(eta)
#' mu.eta <- lf$mu.eta  # this is equivalent to function(eta) exp(eta)
#'
#' fit <- gldrm(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width + Species,
#'              data=iris, linkfun, linkinv, mu.eta)
#'
#' ### Wald 95% confidence interval for Sepal.Width
#' ci <- gldrmCI(fit, "Sepal.Width", test="Wald", level=.95, type="2-sided")
#' ci
#'
#' @export
gldrmCI <- function(gldrmFit, term, test=c("Wald", "LRT"), level=.95,
                    type=c("2-sided", "lb", "ub"), eps=1e-10, maxiter=100)
{
    test <- match.arg(test)
    type <- match.arg(type)
    if (class(gldrmFit) != "gldrmFit")
        stop("gldrmFit should be an object of class gldrmFit, returned from the gldrm function.")
    if (!(term %in% names(gldrmFit$beta)))
        stop("term should be the name of a name from the fitted coefficient vector.")
    if (!(level>0 && level<1))
        stop("level should be between zero and one.")
    if (!(eps > 0))
        stop("eps should be a small value greater than zero.")

    mf <- model.frame(gldrmFit$formula, gldrmFit$data)
    x <- stats::model.matrix(attr(mf, "terms"), mf)
    attributes(x)[c("assign", "contrasts")] <- NULL
    y <- stats::model.response(mf, type="numeric")
    offset <- gldrmFit$object
    linkfun <- gldrmFit$linkfun
    linkinv <- gldrmFit$linkinv
    mu.eta <- gldrmFit$mu.eta
    mu0 <- gldrmFit$mu0
    f0 <- gldrmFit$f0
    offset <- gldrmFit$offset
    beta <- gldrmFit$beta
    seBeta <- gldrmFit$seBeta
    llik <- gldrmFit$llik
    df2 <- gldrmFit$lr.df[2]
    id <- match(term, names(beta))
    cilo <- cihi <- NA  # initialze to NA for one-sided intervals
    iterlo <- iterhi <- NA  # initialize to NA for Wald or one-sided intervals
    pvallo <- pvalhi <- NA  # initialize to NA for Wald or one-sided intervals

    if (type == "2-sided") {
        pvalTarget <- (1 - level) / 2
        step <- stats::qt((level + 1) / 2, df2) * seBeta[id]  # length 2-sided of Wald CI
    } else {
        pvalTarget <- 1 - level
        step <- stats::qt(level, df2) * seBeta[id]  # length of 1-sided Wald CI
    }

    ### Wald test CI
    if (test == "Wald" && type %in% c("2-sided", "lb"))
        cilo <- unname(beta[id]) - step
    if (test == "Wald" && type %in% c("2-sided", "ub"))
        cihi <- unname(beta[id]) + step

    ### Likelihood ratio CI
    if (test=="LRT" && type %in% c("2-sided", "lb")) {
        cilo.hi <- unname(beta[id])
        cilo.lo <- -Inf
        pvallo.hi <- 1
        pvallo.lo <- 0
        betaStart <- beta[-id]
        f0Start <- f0
        stepfact <- 1
        iterlo <- 0
        conv <- FALSE
        while (!conv && iterlo<maxiter) {
            iterlo <- iterlo + 1

            if (cilo.lo == -Inf) {
                cilo.try <- unname(beta[id]) - stepfact * step
                stepfact <- stepfact * 2
            } else {
                # wt <- (pvalTarget - pvallo.lo) / (pvallo.hi - pvallo.lo)
                wt <- .5
                cilo.try <- (1-wt) * cilo.lo + wt * cilo.hi
            }

            ## Use current beta values as starting values unless they violate the convex hull condition.
            fit <- tryCatch({
                gldrm(y ~ x[, -id, drop=FALSE] - 1, data=NULL,
                      linkfun=linkfun, linkinv=linkinv, mu.eta=mu.eta,
                      mu0=mu0, offset=offset + cilo.try * x[, id],
                      gldrmControl=gldrm.control(f0Start=f0Start, betaStart=betaStart))
            }, error=function(e) {
                gldrm(y ~ x[, -id, drop=FALSE] - 1, data=NULL,
                      linkfun=linkfun, linkinv=linkinv, mu.eta=mu.eta,
                      mu0=mu0, offset=offset + cilo.try * x[, id],
                      gldrmControl=gldrm.control(f0Start=f0Start))
            })

            betaStart <- fit$beta
            f0Start <- fit$f0
            pvallo.try <- 1 - stats::pf(2 * (llik - fit$llik), 1, df2)
            conv <- abs(pvallo.try -  pvalTarget) < eps

            if (pvallo.try > pvalTarget) {
                cilo.hi <- cilo.try
                pvallo.hi <- pvallo.try
            } else {
                cilo.lo <- cilo.try
                pvallo.lo <- pvallo.try
            }
        }
        cilo <- cilo.try
        pvallo <- pvallo.try
    }

    if (test=="LRT" && type %in% c("2-sided", "ub")) {
        cihi.hi <- Inf
        cihi.lo <- unname(beta[id])
        pvalhi.hi <- 0
        pvalhi.lo <- 1
        betaStart <- beta[-id]
        f0Start <- f0
        stepfact <- 1
        iterhi <- 0
        conv <- FALSE
        while (!conv && iterhi<maxiter) {
            iterhi <- iterhi + 1

            if (cihi.hi == Inf) {
                cihi.try <- unname(beta[id]) + stepfact * step
                stepfact <- stepfact * 2
            } else {
                # wt <- (pvalTarget - pvalhi.hi) / (pvalhi.lo - pvalhi.hi)
                wt <- .5
                cihi.try <- wt * cihi.lo + (1-wt) * cihi.hi
            }

            ## Use current beta values as starting values unless they violate the convex hull condition.
            fit <- tryCatch({
                gldrm(y ~ x[, -id, drop=FALSE] - 1, data=NULL,
                      linkfun=linkfun, linkinv=linkinv, mu.eta=mu.eta,
                      mu0=mu0, offset=offset + cihi.try * x[, id],
                      gldrmControl=gldrm.control(f0Start=f0Start, betaStart=betaStart))
            }, error=function(e) {
                gldrm(y ~ x[, -id, drop=FALSE] - 1, data=NULL,
                      linkfun=linkfun, linkinv=linkinv, mu.eta=mu.eta,
                      mu0=mu0, offset=offset + cihi.try * x[, id],
                      gldrmControl=gldrm.control(f0Start=f0Start))
            })

            betaStart <- fit$beta
            f0Start <- fit$f0
            pvalhi.try <- 1 - stats::pf(2 * (llik - fit$llik), 1, df2)
            conv <- abs(pvalhi.try -  pvalTarget) < eps

            if (pvalhi.try > pvalTarget) {
                cihi.lo <- cihi.try
                pvalhi.lo <- pvalhi.try
            } else {
                cihi.hi <- cihi.try
                pvalhi.hi <- pvalhi.try
            }
        }
        cihi <- cihi.try
        pvalhi <- pvalhi.try
    }

    ci <- list(term=term, test=test, level=level, type=type, cilo=cilo, cihi=cihi,
               iterlo=iterlo, iterhi=iterhi, pvallo=pvallo, pvalhi=pvalhi)
    class(ci) <- "gldrmCI"
    ci
}

#' Print confidence interval
#'
#' Print method for gldrmCI objects.
#'
#' @param x An S3 object of class 'gldrmCI'.
#' @param digits Number of digits for rounding.
#' @param ... Not used. Additional arguments for print method.
#'
#' @export
print.gldrmCI <- function(x, digits=3, ...)
{
    level <- round(x$level * 100)
    fmt <- paste0("%.", digits, "f")
    cilo <- sprintf(fmt, x$cilo)
    cihi <- sprintf(fmt, x$cihi)

    if (x$test == "Wald") test <- "Wald"
    if (x$test == "LRT") test <- "likelihood ratio"

    if (x$type == "2-sided") {
        cat('\n', level, "% ", test, " confidence interval for ", x$term, ":\n", sep='')
        cat("    (", cilo, ", ", cihi, ")\n", sep='')
    }

    if (x$type == "lb") {
        cat('\n', level, "% ", test, " lower bound for ", x$term, ":\n", sep='')
        cat("    ", cilo, "\n", sep='')
    }

    if (x$type == "ub") {
        cat('\n', level, "% ", test, " upper bound for ", x$term, ":\n", sep='')
        cat("    ", cihi, "\n", sep='')
    }

    return(NULL)
}
