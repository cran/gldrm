% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getTheta.R
\name{theta.control}
\alias{theta.control}
\title{Control arguments for \eqn{\theta} update algorithm}
\usage{
theta.control(
  eps = 1e-10,
  maxiter = 100,
  maxhalf = 20,
  maxtheta = 500,
  logit = TRUE,
  logsumexp = FALSE
)
}
\arguments{
\item{eps}{Convergence threshold for theta updates. Convergence is
evaluated separately for each observation. An observation has converged when
the difference between \eqn{b'(\theta)} and \eqn{\mu} is less than \code{epsTheta}.}

\item{maxiter}{Maximum number of iterations.}

\item{maxhalf}{Maximum number of half steps allowed per iteration if the
convergence criterion does not improve.}

\item{maxtheta}{Absolute value of theta is not allowed to exceed \code{maxtheta}.}

\item{logit}{Logical for whether logit transformation should be used. Use of
this stabilizing transformation appears to be faster in general. Default is TRUE.}

\item{logsumexp}{Logical argument for whether log-sum-exp trick should be used.
This may improve numerical stability at the expense of computational time.}
}
\value{
Object of S3 class "thetaControl", which is a list of control arguments.
}
\description{
This function returns control arguments for the \eqn{\theta} update algorithm.
Each argument has a default value, which will be used unless a different
value is provided by the user.
}
