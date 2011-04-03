\name{K.table}
\title{Tables of K-factors for Tolerance Intervals Based on Normality}
\alias{K.table}
\usage{
K.table(n, alpha, P, side = 1, f = NULL,
        by.arg = c("n", "alpha", "P")) 
}
\description{
  Tabulated summary of k-factors for tolerance intervals based on normality.  The user can specify multiple values
  for each of the three inputs.
}

\arguments{
  \item{n}{A vector of (effective) sample sizes.}
  \item{alpha}{The level chosen such that \code{1-alpha} is the confidence level.  Can be a vector.}
  \item{P}{The proportion of the population to be covered by this tolerance interval.  Can be a vector.}
  \item{side}{Whether a 1-sided or 2-sided tolerance interval is required (determined by \code{side = 1} or \code{side = 2},
  respectively).}
  \item{f}{The number of degrees of freedom associated with calculating the estimate of the population standard deviation.
  If \code{NULL}, then \code{f} is taken to be \code{n-1}.  Only a single value can be specified for \code{f}.}
  \item{by.arg}{How you would like the output organized.  If \code{by.arg = "n"}, then the output provides a list of matrices
  sorted by the values specified in \code{n}.  The matrices have rows corresponding to the values specified by \code{1-alpha} and columns
  corresponding to the values specified by \code{P}.  If \code{by.arg = "alpha"}, then the output provides a list of matrices
  sorted by the values specified in \code{1-alpha}.  The matrices have rows corresponding to the values specified by \code{n} and columns
  corresponding to the values specified by \code{P}.  If \code{by.arg = "P"}, then the output provides a list of matrices
  sorted by the values specified in \code{P}.  The matrices have rows corresponding to the values specified by \code{1-alpha} and columns
  corresponding to the values specified by \code{n}. }
} \value{
  \code{K.table} returns a list with a structure determined by the argument \code{by.arg} described above.
} \seealso{
  \code{\link{K.factor}}
} \details{
  The method used for estimating the k-factors is that due to Howe as it is generally viewed as more accurate than the Weissberg-Beatty method.
}

\references{
  Howe, W. G. (1969), Two-Sided Tolerance Limits for Normal Populations - Some Improvements, \emph{Journal of the
  American Statistical Association}, \bold{64}, 610--620.

  Weissberg, A. and Beatty, G. (1969), Tables of Tolerance Limit Factors for Normal Distributions, \emph{Technometrics},
  \bold{2}, 483--500.
} 

\examples{ 
## Tables generated for each value of the sample size.

K.table(n = seq(50, 100, 10), alpha = c(0.01, 0.05, 0.10), 
        P = c(0.90, 0.95, 0.99), by.arg = "n")

## Tables generated for each value of the confidence level.

K.table(n = seq(50, 100, 10), alpha = c(0.01, 0.05, 0.10), 
        P = c(0.90, 0.95, 0.99), by.arg = "alpha")

## Tables generated for each value of the coverage proportion.

K.table(n = seq(50, 100, 10), alpha = c(0.01, 0.05, 0.10), 
        P = c(0.90, 0.95, 0.99), by.arg = "P")
        
}

\keyword{file}



