\name{uniftol.int}
\title{Uniform Tolerance Intervals}
\alias{uniftol.int}
\usage{
uniftol.int(x, alpha = 0.05, P = 0.99)
}
\description{
  Provides 1-sided tolerance intervals for data distributed according to a uniform distribution.
}

\arguments{
  \item{x}{A vector of data which is distributed according to a uniform distribution.}
  \item{alpha}{The level chosen such that \code{1-alpha} is the confidence level.}
  \item{P}{The proportion of the population to be covered by this tolerance interval.}
} \value{
  \code{uniftol.int} returns a data frame with items:
  \item{alpha}{The specified significance level.}
  \item{P}{The proportion of the population covered by this tolerance interval.}
  \item{1-sided.lower}{The 1-sided lower tolerance bound.}
  \item{1-sided.upper}{The 1-sided upper tolerance bound.}
} 

\references{
  Faulkenberry, G. D. and Weeks, D. L. (1968), Sample Size Determination for Tolerance Limits,
  \emph{Technometrics}, \bold{10}, 343--348.
} 

\examples{ 

## 90\%/90\% 1-sided uniform tolerance intervals for a sample
## of size 50. 

set.seed(100)
x <- runif(50, 0, 50)
out <- uniftol.int(x = x, alpha = 0.10, P = 0.90) 
out

plottol(out, x, plot.type = "hist", side = "two", 
        x.lab = "Uniform Data")
}

\keyword{file}





