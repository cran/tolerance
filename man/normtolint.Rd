\name{normtol.int}
\title{Normal (or Log-Normal) Tolerance Intervals}
\alias{normtol.int}
\usage{
normtol.int(x, alpha = 0.05, P = 0.99, side = 1,
            method = c("HE", "WBE"), log.norm = FALSE)
}
\description{
  Provides 1-sided or 2-sided tolerance intervals for data distributed according to either a normal
  distribution or log-normal distribution.
}

\arguments{
  \item{x}{A vector of data which is distributed according to either a normal distribution or a log-normal distribution.}
  \item{alpha}{The level chosen such that \code{1-alpha} is the confidence level.}
  \item{P}{The proportion of the population to be covered by this tolerance interval.}
  \item{side}{Whether a 1-sided or 2-sided tolerance interval is required (determined by \code{side = 1} or \code{side = 2},
  respectively).}
  \item{method}{The method for calculating the k-factors.  The k-factor for the 1-sided tolerance intervals
  is performed exactly and thus the same for either method chosen.  \code{"HE"} is the
  Howe method and is often viewed as being extremely accurate, even for small sample sizes. \code{"WBE"} is the 
  Weissberg-Beatty method, which performs similarly to the Howe method for larger sample sizes.}
  \item{log.norm}{If \code{TRUE}, then the data is considered to be from a log-normal distribution, in which
  case the output gives tolerance intervals for the log-normal distribution.  The default is \code{FALSE}.} 
} \value{
  \code{normtol.int} returns a data frame with items:
  \item{alpha}{The specified significance level.}
  \item{P}{The proportion of the population covered by this tolerance interval.}
  \item{x.bar}{The sample mean.}
  \item{1-sided.lower}{The 1-sided lower tolerance bound.  This is given only if \code{side = 1}.}
  \item{1-sided.upper}{The 1-sided upper tolerance bound.  This is given only if \code{side = 1}.}
  \item{2-sided.lower}{The 2-sided lower tolerance bound.  This is given only if \code{side = 2}.}
  \item{2-sided.upper}{The 2-sided upper tolerance bound.  This is given only if \code{side = 2}.}
} \seealso{
  \code{\link{Normal}}, \code{\link{K.factor}}
}

\details{
Recall that if the random variable \eqn{X} is distributed according to a log-normal distribution, then the random variable \eqn{Y = ln(X)} is
distributed according to a normal distribution.
}


\references{
  Howe, W. G. (1969), Two-Sided Tolerance Limits for Normal Populations - Some Improvements, \emph{Journal of the
  American Statistical Association}, \bold{64}, 610--620.

  Weissberg, A. and Beatty, G. (1969), Tables of Tolerance Limit Factors for Normal Distributions, \emph{Technometrics},
  \bold{2}, 483--500.
} 

\examples{ 

## 95\%/95\% 2-sided normal tolerance intervals for a sample
## of size 100. 

set.seed(100)
x <- rnorm(100, 0, 0.2)
out <- normtol.int(x = x, alpha = 0.05, P = 0.95, side = 2,
                   method = "HE", log.norm = FALSE)
out

plottol(out, x, plot.type = "both", side = "two", 
        x.lab = "Normal Data")
}

\keyword{file}




