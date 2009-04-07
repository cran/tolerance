\name{K.factor}
\title{Estimating K-factors for Tolerance Intervals Based on Normality}
\alias{K.factor}
\usage{
K.factor(n, f = NULL, alpha = 0.05, P = 0.99, side = 1, 
         method = c("HE", "WBE", "ELL"))
}
\description{
  Estimates k-factors for tolerance intervals based on normality.
}

\arguments{
  \item{n}{The (effective) sample size.}
  \item{f}{The number of degrees of freedom associated with calculating the estimate of the population standard deviation.
  If \code{NULL}, then \code{f} is taken to be \code{n-1}.}
  \item{alpha}{The level chosen such that \code{1-alpha} is the confidence level.}
  \item{P}{The proportion of the population to be covered by the tolerance interval.}
  \item{side}{Whether a 1-sided or 2-sided tolerance interval is required (determined by \code{side = 1} or \code{side = 2},
  respectively).}
  \item{method}{The method for calculating the k-factors.  The k-factor for the 1-sided tolerance intervals
  is performed exactly and thus is the same for the chosen method.  \code{"HE"} is the
  Howe method and is often viewed as being extremely accurate, even for small sample sizes. \code{"WBE"} is the 
  Weissberg-Beatty method, which performs similarly to the Howe method for larger sample sizes. \code{"ELL"} is
  the Ellison correction to the Weissberg-Beatty method when \code{f} is appreciably larger than \code{n^2}. A warning
  message is displayed if \code{f} is not larger than \code{n^2}.}
} \value{
  \code{K.factor} returns the k-factor for tolerance intervals based on normality with the arguments specified above.
} \seealso{
  \code{\link{K.table}}, \code{\link{normtol.int}}
}

\references{
  Ellison, B. E. (1964), On Two-Sided Tolerance Intervals for a Normal Distribution, \emph{Annals of Mathematical
  Statistics}, \bold{35}, 762--772.
  
  Howe, W. G. (1969), Two-Sided Tolerance Limits for Normal Populations - Some Improvements, \emph{Journal of the
  American Statistical Association}, \bold{64}, 610--620.

  Weissberg, A. and Beatty, G. (1969), Tables of Tolerance Limit Factors for Normal Distributions, \emph{Technometrics},
  \bold{2}, 483--500.
} 

\examples{ 

## Showing the effect of the Howe and Weissberg-Beatty 
## estimation methods as the sample size increases.

K.factor(10, P = 0.95, side = 2, method = "HE")
K.factor(10, P = 0.95, side = 2, method = "WBE")

K.factor(100, P = 0.95, side = 2, method = "HE")
K.factor(100, P = 0.95, side = 2, method = "WBE")

K.factor(1000, P = 0.95, side = 2, method = "HE")
K.factor(1000, P = 0.95, side = 2, method = "WBE")

}

\keyword{file}


