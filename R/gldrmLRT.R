#' Likelihood ratio test for nested models
#'
#' Performs a likelihood ratio F-test between nested gldrm models.
#' The F-statistic is calculated as \eqn{2 \times (llik - llik_0) / r}, where
#' \eqn{r} is the difference is the number of parameters between the full and null
#' models. The F-statistic has degrees of freedom \eqn{r} and \eqn{n-p}, where
#' \eqn{n} is the number of observations and \eqn{p} is the number of parameters
#' in the full model.
#'
#' @param gldrmFit The full model. An object of S3 class 'gldrmFit' returned from
#' the gldrm function.
#' @param gldrmFitNull The sub-model being tested under the null hypotheses.
#' An object of S3 class 'gldrmFit' returned from the gldrm function.
#'
#' @return An S3 object of class 'gldrmLRT', containing numerator and denominator
#' degrees of freedom, an F-statistic, and a p-value.
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
#' ### Fit gldrm without the categorical variable "Species"
#' fit0 <- gldrm(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width,
#'               data=iris, linkfun, linkinv, mu.eta)
#'
#' ### Likelihood ratio test for the nested models
#' lrt <- gldrmLRT(fit, fit0)
#' lrt
#'
#' @export
gldrmLRT <- function(gldrmFit, gldrmFitNull)
{
    beta <- gldrmFit$beta
    beta0 <- gldrmFitNull$beta
    p <- sum(!is.na(beta))
    p0 <- sum(!is.na(beta0))
    n <- length(gldrmFit$mu)
    n0 <- length(gldrmFitNull$mu)
    llik <- gldrmFit$llik
    llik0 <- gldrmFitNull$llik

    if (class(gldrmFit) != "gldrmFit")
        stop("gldrmFit must be an S3 object of class gldrmFit.")
    if (class(gldrmFitNull) != "gldrmFit")
        stop("gldrmFitNull must be an S3 object of class gldrmFit.")
    if (p0 >= p)
        stop("gldrmFitNull must be a sub-model of gldrmFit")
    if (n != n0)
        stop("gldrm and gldrmNull have a different number of observations.")
    if (!all(names(beta0) %in% names(beta)))
        warning(paste0("Coefficient names of the nested model are not a subset of names ",
                       "in the full model. Make sure the models are nested."))

    df <- c(p-p0, n-p)
    fstat <- 2 * (llik - llik0) / df[1]
    pval <- 1 - stats::pf(fstat, df[1], df[2])

    lrt <- list(df=df, fstat=fstat, pval=pval)
    class(lrt) <- "gldrmLRT"
    lrt
}

#' Print likelihood ratio test results
#'
#' Print method for gldrmLRT objects. Prints results of a likelihood ratio F-test
#' between nested models.
#'
#' @param x An S3 object of class 'gldrmLRT'.
#' @param digits Number of digits for rounding.
#' @param ... Not used. Additional arguments for print method.
#'
#' @export
print.gldrmLRT <- function(x, digits=3, ...)
{
    if (x$pval < 2e-16) {
        pval <- "< 2e-16"
    } else {
        pval <- signif(x$pval, digits)
    }

    cat("\nLikelihood ratio test:\n\n")
    cat("                   F-statistic: ", signif(x$fstat, digits), "\n")
    cat("  Numerator degrees of freedom: ", x$df[1], "\n")
    cat("Denomicator degrees of freedom: ", x$df[2], "\n")
    cat("                       P-value: ", pval, "\n")

    return(NULL)
}
