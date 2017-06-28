#' Control arguments for \eqn{\beta} update algorithm
#'
#' This function returns control arguments for the \eqn{\beta} update algorithm.
#' Each argument has a default value, which will be used unless a different
#' value is provided by the user.
#'
#' @param eps Convergence threshold. The update has converged when the relative
#' change in log-likelihood between iterations is less than \code{eps}.
#' Only applies if \code{maxiter>1}.
#' @param maxiter Maximum number of iterations allowed.
#' @param maxhalf Maximum number of half steps allowed per iteration if
#' log-likelihood does not improve.
#'
#' @return Object of S3 class "betaControl", which is a list of control arguments.
#'
#' @export
beta.control <- function(eps=1e-10, maxiter=1, maxhalf=10)
{
    betaControl <- as.list(environment())
    class(betaControl) <- "betaControl"
    betaControl
}

#' @keywords internal
getBeta <- function(betaStart, thStart, llik, f0,
                    spt, y, ySptIndex, x, linkinv, mu.eta, offset,
                    thetaControl=theta.control(), betaControl=beta.control())
{
    ## Extract control arguments
    if (class(betaControl) != "betaControl")
      stop("betaControl must be an object of class betaControl returned by betaControl() function.")
    eps <- betaControl$eps
    maxiter <- betaControl$maxiter
    maxhalf <- betaControl$maxhalf

    sptMin <- min(spt)
    sptMax <- max(spt)
    beta <- betaStart
    th <- thStart

    conv <- FALSE
    maxhalfreached <- FALSE
    iter <- 0
    while (!conv && !maxhalfreached && iter < maxiter)
    {
        iter <- iter+1

        ## Update mean vector and related quantities
        eta <- c(x %*% beta + offset)
        mu <- linkinv(eta)
        dmudeta <- mu.eta(eta)
        betaold <- beta
        muold <- mu
        thold <- th
        llikold <- llik

        ## Compute weighted least squares update
        varY <- th$bPrime2
        wSqrt <- dmudeta / sqrt(varY)
        wSqrt[dmudeta==0] <- 0
        if (any(wSqrt==Inf)) break
        # z <- (y-mu) / dmudeta + eta (don't add eta: want beta step, not new beta)
        r <- (y-mu) / dmudeta
        # betastep <- stats::lm.wfit(x, r, wSqrt^2)$coef
        betastep <- qr.coef(qr(wSqrt*x), wSqrt*r)
        betastep[is.na(betastep)] <- 0

        ### Update beta and take half steps if log-likelihood does not improve
        beta <- beta + betastep
        eta <- c(x %*% beta + offset)
        mu <- linkinv(eta)
        if (min(mu)<sptMin || max(mu)>sptMax) {
            llik <- -Inf
        } else {
            th <- getTheta(spt=spt, f0=f0, mu=mu, thetaStart=thold$theta, thetaControl=thetaControl)
            llik <- getLoglik(ySptIndex=ySptIndex, fTilt=th$fTilt)
        }

        nhalf <- 0
        while ((llik<llikold) && (nhalf<maxhalf)) {
            nhalf <- nhalf + 1
            beta <- (beta + betaold) / 2
            eta <- c(x %*% beta + offset)
            mu <- linkinv(eta)
            if (min(mu)<sptMin || max(mu)>sptMax) {
                llik <- -Inf
            } else {
                th <- getTheta(spt=spt, f0=f0, mu=mu, thetaStart=thold$theta, thetaControl=thetaControl)
                llik <- getLoglik(ySptIndex=ySptIndex, fTilt=th$fTilt)
            }
        }

        if (llik < llikold) {
            beta <- betaold
            mu <- muold
            th <- thold
            llik <- llikold
            conv <- FALSE
            maxhalfreached <- TRUE
        } else {
            del <- (llik - llikold) / llik
            if (llik == 0) del <- 0  # consider converged if model fit is perfect
            conv <- del < eps
        }
    }

    return(list(beta=beta, mu=mu, th=th, llik=llik, iter=iter, conv=conv))
}
