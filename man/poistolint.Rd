\name{poistol.int}
\title{Poisson Tolerance Intervals}
\alias{poistol.int}
\usage{
poistol.int(x, n, m, alpha = 0.05, P = 0.99, side = 1, 
            method = c("TAB", "LS", "SC"))
}
\description{
  Provides 1-sided or 2-sided tolerance intervals for Poisson random variables.  From a statistical quality control
  perspective, these limits bound the number of occurrences (which follow a Poisson distribution) in a specified future time period.
}

\arguments{
  \item{x}{The number of occurrences of the event in time period \code{n}. Can be a vector of length \code{n}, in which case the sum of \code{x} is used.}
  \item{n}{The time period of the original measurements.}
  \item{m}{The specified future length of time.}
  \item{alpha}{The level chosen such that \code{1-alpha} is the confidence level.}
  \item{P}{The proportion of occurrences in future time lengths of size \code{m} 
  to be covered by this tolerance interval.}
  \item{side}{Whether a 1-sided or 2-sided tolerance interval is required (determined by \code{side = 1} or \code{side = 2},
  respectively).}
  \item{method}{The method for calculating the lower and upper confidence bounds, which are used in the calculation
  of the tolerance bounds.  The default method is \code{"TAB"}, which is the tabular method and is usually preferred for a smaller
  number of occurrences.  \code{"LS"} gives the large-sample method, which is usually preferred when the number of occurrences is
  \code{x>20}. code{"SC"} gives the score method, which again is usually used when the number of occurrences is relatively large.}

} \value{
  \code{poistol.int} returns a data frame with items:
  \item{alpha}{The specified significance level.}
  \item{P}{The proportion of occurrences in future time periods of length \code{m}.}
  \item{lambda.hat}{The mean occurrence rate per unit time, calculated by \code{x/n}.}
  \item{1-sided.lower}{The 1-sided lower tolerance bound.  This is given only if \code{side = 1}.}
  \item{1-sided.upper}{The 1-sided upper tolerance bound.  This is given only if \code{side = 1}.}
  \item{2-sided.lower}{The 2-sided lower tolerance bound.  This is given only if \code{side = 2}.}
  \item{2-sided.upper}{The 2-sided upper tolerance bound.  This is given only if \code{side = 2}.}
} \seealso{
  \code{\link{Poisson}}, \code{\link{umatol.int}}
}

\references{
  Hahn, G. J. and Chandra, R. (1981), Tolerance Intervals for Poisson and Binomial Variables,
  \emph{Journal of Quality Technology}, \bold{13}, 100--110.
} 

\examples{ 
## 95%/90% 1-sided Poisson tolerance limits for future 
## occurrences in a period of length 3.  All three methods
## are presented for comparison.

poistol.int(x = 45, n = 9, m = 3, alpha = 0.05, P = 0.90,
            side = 1, method = "TAB")
poistol.int(x = 45, n = 9, m = 3, alpha = 0.05, P = 0.90,
            side = 1, method = "LS")
poistol.int(x = 45, n = 9, m = 3, alpha = 0.05, P = 0.90,
            side = 1, method = "SC")

## 95%/90% 2-sided Poisson tolerance intervals for future 
## occurrences in a period of length 15.  All three methods
## are presented for comparison.

poistol.int(x = 45, n = 9, m = 15, alpha = 0.05, P = 0.90,
            side = 2, method = "TAB")
poistol.int(x = 45, n = 9, m = 15, alpha = 0.05, P = 0.90,
            side = 2, method = "LS")
poistol.int(x = 45, n = 9, m = 15, alpha = 0.05, P = 0.90,
            side = 2, method = "SC")
}

\keyword{file}


