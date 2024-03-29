% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gldrmPIT.R
\name{gldrmPIT}
\alias{gldrmPIT}
\title{Confidence intervals for gldrm coefficients}
\usage{
gldrmPIT(
  gldrmFit,
  nbreaks = 7,
  cex.main = NULL,
  cex.lab = NULL,
  cex.axis = NULL
)
}
\arguments{
\item{gldrmFit}{A gldrm model fit. Must be an S3 object of class "gldrm",
returned from the \code{gldrm} function. The matrix of semiparametric
tilted probabilities must be returned, which is done by fitting gldrm with
\code{gldrmControl = gldrm.control(returnfTiltMatrix = TRUE)}.}

\item{nbreaks}{Number of breaks in the histogram.}

\item{cex.main}{Text size for main titles.}

\item{cex.lab}{Text size for axis labels.}

\item{cex.axis}{Text size for axis numbers.}
}
\value{
Randomized robability inverse transform as a vector. Also plots the 
histogram and uniform QQ plot.
}
\description{
Plots and returns the randomized probability inverse transform of a 
fitted gldrm.
}
\details{
The probability inverse transform is defined generally as \eqn{\hat{F}(y|x)},
which is the fitted conditional cdf of each observation evaluated at the
observed response value. In the case of gldrm, the fitted cdf is descrete, so
we draw a random value from a uniform distribution on the interval 
(\eqn{\hat{F}(y|x)}, \eqn{\hat{F}(y-|x)}), where \eqn{y-} is the next largest 
observed support less than \eqn{y} (or -Infinity if \eqn{y} is the minimum 
support value). The output and plots generated by this function will vary
slightly each time it is called (unless the random number generator seed is 
set beforehand).
}
\examples{
data(iris, package="datasets")

### Fit gldrm and return fTiltMatrix
fit <- gldrm(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width + Species,
             data=iris, link="log")
             
# Probability inverse transform plot
gldrmPIT(fit)

}
