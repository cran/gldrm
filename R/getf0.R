################################################################################
# These functions are only used for the 2 term approximate information
# (method = "approx2")
################################################################################

## Computes inverse of a 2x2 matrix
invert2by2 <- function(m) {
    matrix(c(m[4], -m[2], -m[3], m[1]), nrow=2) /
        (m[1] * m[4] - m[2] * m[3])
}

## Computes (A+BCB')^{-1}, where Ainv is available and B is rank 1 or 2
## Adiag is an indicator of whether A is a diagonal matrix
woodbury <- function(Ainv, Cinv, B) {
    AinvB <- Ainv %*% B
    mid <- Cinv + crossprod(B, AinvB)
    midinv <- invert2by2(mid)
    inv <- Ainv - AinvB %*% tcrossprod(midinv, AinvB)
    inv
}
################################################################################

#' Control arguments for f0 update algorithm
#'
#' This function returns control arguments for the \eqn{f_0} update algorithm.
#' Each argument has a default value, which will be used unless a different
#' value is provided by the user.
#'
#' @param eps Convergence threshold. The update has converged when the relative
#' change in log-likelihood between iterations is less than \code{eps}.
#' absolute change is less than \code{thesh}.
#' @param maxiter Maximum number of iterations allowed.
#' @param maxhalf Maximum number of half steps allowed per iteration if
#' log-likelihood does not improve between iterations.
#' @param maxlogstep Maximum optimization step size allowed on the
#' \code{log(f0)} scale.
#'
#' @return Object of S3 class "f0Control", which is a list of control arguments.
#'
#' @export
f0.control <- function(eps=1e-10, maxiter=1000, maxhalf=20, maxlogstep=2)
{
    f0Control <- as.list(environment())
    class(f0Control) <- "f0Control"
    f0Control
}

#' @keywords internal
getf0 <- function(f0Start, thStart, llik, mu, mu0, ySptIndex, spt,
                  thetaControl=theta.control(), f0Control=f0.control(),
                  method=c("bfgs", "approx2", "gradient"), trace=FALSE)
{
    method <- match.arg(method)  # gldrm function uses bfgs, the default

    ## Extract theta control arguments
    if (class(f0Control) != "f0Control")
        stop("f0Control must be an object of class f0Control returned by f0Control() function.")
    eps <- f0Control$eps
    maxiter <- f0Control$maxiter
    maxhalf <- f0Control$maxhalf
    maxlogstep <- f0Control$maxlogstep

    f0 <- f0Start
    f0 <- f0 / sum(f0)
    f0 <- getTheta(spt=spt, f0=f0, mu=mu0, thetaStart=0, thetaControl=thetaControl)$fTilt[, 1]
    th <- thStart
    y <- spt[ySptIndex]
    yFreq <- table(ySptIndex)
    attributes(yFreq) <- NULL
    ymm <- outer(spt, mu, "-")
    yeqmu <- which(y == mu)
    score.log <- NULL

    conv <- FALSE
    iter <- 0
    while (!conv && iter<maxiter) {
        iter <- iter + 1

        # Score calculation
        score.logOld <- score.log
        fTiltSums <- rowSums(th$fTilt)
        ymmfTilt <- ymm * th$fTilt
        ystd <- (y - mu) / th$bPrime2
        ystd[yeqmu] <- 0  # prevent 0/0
        score.logT1 <- yFreq
        score.logT2 <- fTiltSums
        score.logT3 <- c(ymmfTilt %*% ystd)
        score.log <- score.logT1 - score.logT2 - score.logT3

        if (method == "bfgs")
        {
            # Inverse info, score step, and f0 step are on the log scale (score is not)
            if (iter == 1) {
                d1 <- min(fTiltSums)  # max inverse diagonal of first information term, on log scale
                d2 <- max(abs(score.log)) / maxlogstep
                d <- max(d1, d2)
                infoinvBFGS.log <- diag(1/d, nrow=length(f0))
            } else {
                scorestep.log <- score.log - score.logOld
                f0step.log <- log(f0) - log(f0old)
                sy <- sum(f0step.log * scorestep.log)
                yiy <- c(crossprod(scorestep.log, infoinvBFGS.log %*% scorestep.log))
                iys <- tcrossprod(infoinvBFGS.log %*% scorestep.log, f0step.log)
                infoinvBFGS.log <- infoinvBFGS.log + ((yiy - sy) / sy^2) * tcrossprod(f0step.log) - (1 / sy) * (iys + t(iys))
            }
            logstep <- c(infoinvBFGS.log %*% score.log)
        }

        if (method == "approx2")
        {
            n <- length(mu)
            fTiltMeans <- fTiltSums / n
            fTiltSD <- sqrt(rowMeans((th$fTilt - fTiltMeans)^2))
            d <- fTiltSums
            B <- cbind(fTiltMeans, fTiltSD)
            Ainv <- diag(1/d)
            Cinv <- diag(rep(-1/n, 2))
            # info <- diag(d) - n * tcrossprod(B)
            infoinv <- woodbury(Ainv, Cinv, B)  # on log scale
            logstep <- c(infoinv %*% score.log)
            infoinvBFGS.log <- NA
        }

        if (method == "gradient") {
            # Save values from previous iteration
            f0old <- f0
            thold <- th
            llikold <- llik

            d1 <- min(fTiltSums)  # max inverse diagonal of first information term, on log scale
            d2 <- max(abs(score.log)) / maxlogstep
            d <- max(d1, d2)
            logstep <- score.log / d
            infoinvBFGS.log <- NA
        }

        # Cap log(f0) step size
        logstep.max <- max(abs(logstep))
        if (logstep.max > maxlogstep)
            logstep <- logstep * (maxlogstep / logstep.max)

        # Save values from previous iteration
        f0old <- f0
        thold <- th
        llikold <- llik

        # Take update step
        f0 <- exp(log(f0) + logstep)
        # Scale and tilt f0
        f0 <- f0 / sum(f0)
        f0 <- getTheta(spt=spt, f0=f0, mu=mu0, thetaStart=0, thetaControl=thetaControl)$fTilt[, 1]
        # Update theta and likelihood
        thold <- th
        llikold <- llik
        th <- getTheta(spt, f0, mu, thetaStart=th$theta, thetaControl=thetaControl)
        llik <- getLoglik(ySptIndex, th$fTilt)
        conv <- abs((llik - llikold) / (llik + 1e-100)) < eps

        # If log-likelihood does not improve, change step direction to be along gradient
        # Take half steps until likelihood improves
        # Continue taking half steps until log likelihood no longer improves
        nhalf <- 0
        if (llik<llikold) {
            llikprev <- -Inf
            while ((llik<llikold || llik>llikprev) && nhalf<maxhalf) {
                nhalf <- nhalf + 1

                # Set previous values
                llikprev <- llik
                thprev <- th
                f0prev <- f0
                infoinvBFGS.logprev <- infoinvBFGS.log

                f0 <- exp((log(f0) + log(f0old)) / 2)
                f0 <- f0 / sum(f0)
                f0 <- getTheta(spt=spt, f0=f0, mu=mu0, thetaStart=0, thetaControl=thetaControl)$fTilt[, 1]
                th <- getTheta(spt, f0, mu, thetaStart=th$theta, thetaControl=thetaControl)
                llik <- getLoglik(ySptIndex, th$fTilt)
                infoinvBFGS.log <- infoinvBFGS.log / 2
            }

            if (llik < llikprev) {
                nhalf <- nhalf - 1
                llik <- llikprev
                th <- thprev
                f0 <- f0prev
                infoinvBFGS.log <- infoinvBFGS.logprev
            }

            conv <- abs((llik - llikold) / (llik + 1e-100)) < eps
        }

        if (llik < llikold) {
            f0 <- f0old
            th <- thold
            llik <- llikold
            conv <- TRUE
        }

        if (trace) {
            printout <- paste0("iter ", iter, ": llik=", llik)
            if (nhalf > 0)
                printout <- paste0(printout, "; ", nhalf, " half steps")
            cat(printout, "\n")
        }
    }

    # Final score calculation
    fTiltSums <- rowSums(th$fTilt)
    ymmfTilt <- ymm * th$fTilt
    ystd <- (y - mu) / th$bPrime2
    ystd[yeqmu] <- 0  # prevent 0/0
    score.logT1 <- yFreq
    score.logT2 <- fTiltSums
    score.logT3 <- c(ymmfTilt %*% ystd)
    score.log <- score.logT1 - score.logT2 - score.logT3

    # Final info calculation
    info.logT1 <- diag(fTiltSums)
    info.logT2 <- tcrossprod(th$fTilt)
    info.logT3 <- tcrossprod(ymmfTilt, ymmfTilt * rep(ystd, each=nrow(ymmfTilt)))
    info.log <- info.logT1 - info.logT2 - info.logT3

    list(f0=f0, llik=llik, th=th, conv=conv, iter=iter, nhalf=nhalf,
         score.log=score.log, info.log=info.log)
}
