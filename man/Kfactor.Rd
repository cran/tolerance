\name{K.factor}
\title{Estimating K-factors for Tolerance Intervals Based on Normality}
\alias{K.factor}
\usage{
K.factor(n, alpha = 0.05, P = 0.99, side = 1, 
         method = c("HE", "WBE"))
}
\description{
  Estimates k-factors for tolerance intervals based on normality by using either the Howe method or the
  Weissberg-Beatty method.
}

\arguments{
  \item{n}{The sample size.}
  \item{alpha}{The level chosen such that \code{1-alpha} is the confidence level.}
  \item{P}{The proportion of the population to be covered by this tolerance interval.}
  \item{side}{Whether a 1-sided or 2-sided tolerance interval is required (determined by \code{side = 1} or \code{side = 2},
  respectively).}
  \item{method}{The method for calculating the k-factors.  The k-factor for the 1-sided tolerance intervals
  is performed exactly and thus the same for either method chosen.  \code{"HE"} is the
  Howe method and is often viewed as being extremely accurate, even for small sample sizes. \code{"WBE"} is the 
  Weissberg-Beatty method, which performs similarly to the Howe method for larger sample sizes.}
} \value{
  \code{K.factor} returns the k-factor for tolerance intervals based on normality with the arguments specified above.
} \seealso{
  \code{\link{K.table}}, \code{\link{normtol.int}}
}

\references{
  Howe, W. G. (1969), Two-Sided Tolerance Limits for Normal Populations - Some Improvements, \emph{Journal of the
  American Statistical Association}, \bold{64}, 610--620.

  Weissberg, A. and Beatty, G. (1969), Tables of Tolerance Limit Factors for Normal Distributions, \emph{Technometrics},
  \bold{2}, 483--500.
} 

\examples{ 

## Showing the effect of the two estimation methods as the
## sample size increases.

K.factor(10, P = 0.95, side = 2, method = "HE")
K.factor(10, P = 0.95, side = 2, method = "WBE")

K.factor(100, P = 0.95, side = 2, method = "HE")
K.factor(100, P = 0.95, side = 2, method = "WBE")

K.factor(1000, P = 0.95, side = 2, method = "HE")
K.factor(1000, P = 0.95, side = 2, method = "WBE")

}

\keyword{file}


