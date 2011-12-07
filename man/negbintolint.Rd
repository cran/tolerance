\name{negbintol.int}
\title{Negative Binomial Tolerance Intervals}
\alias{negbintol.int}
\usage{
negbintol.int(x, n = NULL, N, m, alpha = 0.05, P = 0.99, 
              side = 1)
}
\description{
Provides 1-sided or 2-sided tolerance intervals for negative binomial random variables. From a statistical quality control perspective, these limits use the number of failures that occur to reach \code{N} successes to bound the number of failures for a specified amount of future successes (\code{m}). 
}

\arguments{
  \item{x}{The total number of failures that occur from a sample of size \code{n}.  Can be a vector of length \code{n}, in which case the sum of \code{x} is computed.}
  \item{n}{The size of the random sample from this process.  Note that if the sum of \code{x} is provided, then \code{n} needs to be specified in order to provide correct estimates.}
  \item{N}{The target number of successes (sometimes called size) for each trial.}
  \item{m}{The target number of successes in a future lot for which the tolerance limits will be calculated. If \code{m = NULL}, then the tolerance limits will be constructed assuming \code{N} for the target number of future successes.}
  \item{alpha}{The level chosen such that \code{1-alpha} is the confidence level.}
  \item{P}{The proportion of the defective (or acceptable) units in future samples of size \code{m} 
  to be covered by this tolerance interval.}
  \item{side}{Whether a 1-sided or 2-sided tolerance interval is required (determined by \code{side = 1} or \code{side = 2},
  respectively).}
} \value{
  \code{negbintol.int} returns a data frame with items:
  \item{alpha}{The specified significance level.}
  \item{P}{The proportion of defective (or acceptable) units in future samples of size \code{m}.}
  \item{nu.hat}{The probability of success in each trial, calculated by \code{N/(N+x/n)}.}
  \item{1-sided.lower}{The 1-sided lower tolerance bound.  This is given only if \code{side = 1}.}
  \item{1-sided.upper}{The 1-sided upper tolerance bound.  This is given only if \code{side = 1}.}
  \item{2-sided.lower}{The 2-sided lower tolerance bound.  This is given only if \code{side = 2}.}
  \item{2-sided.upper}{The 2-sided upper tolerance bound.  This is given only if \code{side = 2}.}
} \seealso{
  \code{\link{NegBinomial}}, \code{\link{umatol.int}}
}

\references{
  Hahn, G. J. and Chandra, R. (1981), Tolerance Intervals for Poisson and Binomial Variables,
  \emph{Journal of Quality Technology}, \bold{13}, 100--110.
} 
\details{
This function takes the approach for Poisson and binomial random variables developed in Hahn and Chandra (1981) and applies it to the negative binomial case.
}
\examples{ 
## Comparison of 95%/99% 1-sided tolerance limits with the 
## uniformly most accurate approach when there are 10 trials 
## (resulting in 50 failures) before 3 successes are reached for 
## each trial.

umatol.int(x = 50, n = 10, N = 3, dist = "NegBin")
negbintol.int(x = 50, n = 10, N = 3, side = 1)

## 95%/99% 1-sided tolerance limits and 2-sided tolerance 
## interval for the same setting above, but when we are 
## interested in a future experiment that requires 20 successes 
## be reached for each trial.

negbintol.int(x = 50, n = 10, N = 3, m = 20, side = 1)
negbintol.int(x = 50, n = 10, N = 3, m = 20, side = 2)
}

\keyword{file}

