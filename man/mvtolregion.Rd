\name{mvtol.region}
\title{Multivariate Normal Tolerance Regions}
\alias{mvtol.region}
\usage{
mvtol.region(x, alpha = 0.05, P = 0.99, B = 1000)
}
\description{
  Determines the appropriate tolerance factor for computing multivariate normal tolerance regions based on Monte Carlo simulation.
}

\arguments{
  \item{x}{An \code{n}x\code{p} matrix of data assumed to be drawn from a \code{p}-dimensional multivariate normal
  distribution.  \code{n} pertains to the sample size.}
  \item{alpha}{The level chosen such that \code{1-alpha} is the confidence level.  A vector of \code{alpha} values
  may be specified.}
  \item{P}{The proportion of the population to be covered by this tolerance region.  A vector \code{P} values may be
  specified.}
  \item{B}{The number of iterations used for the Monte Carlo algorithm which determines the tolerance factor.  The number of
  iterations should be at least as large as the default value of 1000.}
} \value{
  \code{mvtol.region} returns a matrix where the rows pertain to each confidence level \code{1-alpha} specified and the columns
  pertain to each proportion level \code{P} specified.
}

\references{
  Krishnamoorthy, K. and Mondal, S. (2006), Improved Tolerance Factors for Multivariate Normal Distributions, 
  \emph{Communications in Statistics - Simulation and Computation}, \bold{35}, 461--478.
} 

\details{
  A basic sketch of how the algorithm works is as follows:
  
  (1) Generate independent chi-square random variables and Wishart random matrices.
  
  (2) Compute the eigenvalues of the randomly generated Wishart matrices.
  
  (3) Iterate the above steps to generate a set of \code{B} sample values such that the \code{100(1-alpha)}-th percentile is an approximate tolerance factor.
}


\examples{ 

## 90\%/90\% bivariate normal tolerance region. 

set.seed(100)
x1 <- rnorm(100, 0, 0.2)
x2 <- rnorm(100, 0, 0.5)
x <- cbind(x1, x2)

out1 <- mvtol.region(x = x, alpha = 0.10, P = 0.90, B = 1000)
out1
plottol(out1, x)

## 90\%/90\% trivariate normal tolerance region. 

set.seed(100)
x1 <- rnorm(100, 0, 0.2)
x2 <- rnorm(100, 0, 0.5)
x3 <- rnorm(100, 5, 1)
x <- cbind(x1, x2, x3)
mvtol.region(x = x, alpha = c(0.10, 0.05, 0.01), 
             P = c(0.90, 0.95, 0.99), B = 1000) 

out2 <- mvtol.region(x = x, alpha = 0.10, P = 0.90, B = 1000) 
out2
plottol(out2, x)
}

\keyword{file}




