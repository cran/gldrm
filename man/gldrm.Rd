% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gldrm.R
\name{gldrm}
\alias{gldrm}
\title{Fits a generalized linear density ratio model (GLDRM)}
\usage{
gldrm(
  formula,
  data = NULL,
  link = "identity",
  mu0 = NULL,
  offset = NULL,
  gldrmControl = gldrm.control(),
  thetaControl = theta.control(),
  betaControl = beta.control(),
  f0Control = f0.control()
)
}
\arguments{
\item{formula}{An object of class "formula".}

\item{data}{An optional data frame containing the variables in the model.}

\item{link}{Link function. Can be a character string to be passed to the
\code{make.link} function in the \code{stats} package (e.g. "identity",
"logit", or "log").
Alternatively, \code{link} can be a list containing three functions named
\code{linkfun}, \code{linkinv}, and \code{mu.eta}. The first is the link
function. The second is the inverse link function. The third is the derivative
of the inverse link function. All three functions must be vectorized.}

\item{mu0}{Mean of the reference distribution. The reference distribution is
not unique unless its mean is restricted to a specific value. This value can
be any number within the range of observed values, but values near the boundary
may cause numerical instability. This is an optional argument with \code{mean(y)}
being the default value.}

\item{offset}{Known component of the linear term. Offset must be passed through
this argument - offset terms in the formula will be ignored.
value and covariate values. If sampling weights are a function of both the
response value and covariates, then \code{sampprobs} must be a \eqn{n \times q}
matrix, where \eqn{n} is the number of observations and \eqn{q} is the number
of unique observed values in the response vector. If sampling weights do not
depend on the covariate values, then \code{sampprobs} may alternatively be passed
as a vector of length \eqn{n}. All values must be nonnegative and are assumed to
correspond to the sorted response values in increasing order.}

\item{gldrmControl}{Optional control arguments.
Passed as an object of class "gldrmControl", which is constructed by the
\code{gldrm.control} function.
See \code{gldrm.control} documentation for details.}

\item{thetaControl}{Optional control arguments for the theta update procedure.
Passed as an object of class "thetaControl", which is constructed by the
\code{theta.control} function.
See \code{theta.control} documentation for details.}

\item{betaControl}{Optional control arguments for the beta update procedure.
Passed as an object of class "betaControl", which is constructed by the
\code{beta.control} function.
See \code{beta.control} documentation for details.}

\item{f0Control}{Optional control arguments for the \code{f0} update procedure.
Passed as an object of class "f0Control", which is constructed by the
\code{f0.control} function.
See \code{f0.control} documentation for details.}
}
\value{
An S3 object of class "gldrm". See details.
}
\description{
A GLDRM is a semiparametric generalized linear model.
In contrast to a GLM, which assumes a particular exponential family distribution,
the GLDRM uses a semiparametric likelihood to estimate the reference distribution.
The reference distribution may be any discrete, continuous, or mixed exponential
family distribution. The model parameters, which include both the regression
coefficients and the cdf of the unspecified reference distribution, are estimated
by maximizing a semiparametric likelihood. Regression coefficients are estimated
with no loss of efficiency, i.e. the asymptotic variance is the same as if the
true exponential family distribution were known.
}
\details{
The arguments \code{linkfun}, \code{linkinv}, and \code{mu.eta}
mirror the "link-glm" class. Objects of this class can be created with the
\code{stats::make.link} function.

The "gldrm" class is a list of the following items.
\itemize{
\item \code{conv} Logical indicator for whether the gldrm algorithm
converged within the iteration limit.
\item \code{iter} Number of iterations used. A single iteration is a \code{beta}
update, followed by an \code{f0} update.
\item \code{llik} Semiparametric log-likelihood of the fitted model.
\item \code{beta} Vector containing the regression coefficient estimates.
\item \code{mu} Vector containing the estimated mean response value for each
observation in the training data.
\item \code{eta} Vector containing the estimated linear combination of
covariates for each observation.
\item \code{f0} Vector containing the semiparametric estimate of the reference
distribution, evaluated at the observed response values. The values of correspond
to the support values, sorted in increasing order.
\item \code{spt} Vector containing the unique observed response values, sorted in
increasing order.
\item \code{mu0} Mean of the estimated semiparametric reference distribution.
The mean of the reference distribution must be fixed at a value in order for
the model to be identifiable. It can be fixed at any value within the range
of observed response values, but the \code{gldrm} function assigns \code{mu0}
to be the mean of the observed response values.
\item \code{varbeta} Estimated variance matrix of the regression coefficients.
\item \code{seBeta} Standard errors for \eqn{\hat{\beta}}. Equal to
\code{sqrt(diag(varbeta))}.
\item \code{seMu} Standard errors for \eqn{\hat{\mu}} computed from \code{varbeta}.
\item \code{seEta} Standard errors for \eqn{\hat{\eta}} computed from \code{varbeta}.
\item \code{theta} Vector containing the estimated tilt parameter for each observation.
The tilted density function of the response variable is given by
\deqn{f(y|x_i) = f_0(y) \exp(\theta_i y) / \int f_0(u) \exp(\theta_i u) du.}
\item \code{bPrime} is a vector containing the mean of the tilted distribution,
\eqn{b'(\theta_i)}, for each observation. \code{bPrime} should match \code{mu},
except in cases where \code{theta} is capped for numerical stability.
\deqn{b'(\theta_i) = \int u f(u|x_i) du}
\item \code{bPrime2} is a vector containing the variance of the tilted
distribution, \eqn{b''(\theta_i)}, for each observation.
\deqn{b''(\theta_i) = \int (u - b'(\theta_i))^2 f(u|x_i) du}
\item \code{fTilt} is a vector containing the semiparametric fitted probability,
\eqn{\hat{f}(y_i | x_i)}, for each observation. The semiparametric
log-likelihood is equal to
\deqn{\sum_{i=1}^n \log \hat{f}(y_i | x_i).}
\item \code{sampprobs} If sampling probabilities were passed through the
\code{sampprobs} argument, then they are returned here in matrix form.
Each row corresponds to an observation.
\item \code{llikNull} Log-likelihood of the null model with no covariates.
\item \code{lr.stat} Likelihood ratio test statistic comparing fitted model to
the null model. It is calculated as \eqn{2 \times (llik - llik_0) / (p-1)}.
The asymptotic distribution is F(p-1, n-p) under the null hypothesis.
\item \code{lr.pval} P-value of the likelihood ratio statistic.
\item \code{fTiltMatrix} is a matrix containing the semiparametric density for
each observation, i.e. \eqn{\hat{f}(y | x_i)} for each unique \code{y} value.
This is a matrix with nrow equal to the number of observations and ncol equal
to the number of unique response values observed.
Only returned if \code{returnfTilt = TRUE} in the gldrmControl arguments.
\item \code{score.logf0} Score function for \code{log(f0)}.
Only returned if \code{returnf0ScoreInfo = TRUE} in the gldrmControl arguments.
\item \code{info.logf0} Information matrix for \code{log(f0)}.
Only returned if \code{returnf0ScoreInfo = TRUE} in the gldrmControl arguments.
\item \code{formula} Model formula.
\item \code{data} Model data frame.
\item \code{link} Link function. If a character string was passed to the
\code{link} argument, then this will be an object of class "link-glm".
Otherwise, it will be the list of three functions passed to the \code{link} argument.
}
}
\examples{
data(iris, package="datasets")

# Fit a gldrm with log link
fit <- gldrm(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width + Species,
             data=iris, link="log")
fit

# Fit a gldrm with custom link function
link <- list()
link$linkfun <- function(mu) log(mu)^3
link$linkinv <- function(eta) exp(eta^(1/3))
link$mu.eta <- function(eta) exp(eta^(1/3)) * 1/3 * eta^(-2/3)
fit2 <- gldrm(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width + Species,
              data=iris, link=link)
fit2

}
