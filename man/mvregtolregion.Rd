\name{mvregtol.region}
\title{Multivariate (Multiple) Linear Regression Tolerance Regions}
\alias{mvregtol.region}
\usage{
mvregtol.region(y, x, new.x = NULL, int = TRUE, alpha = 0.05, 
                P = 0.99, B = 1000)
}
\description{
  Determines the appropriate tolerance factor for computing multivariate (multiple) linear regression tolerance regions based on Monte Carlo simulation.
}

\arguments{
  \item{y}{An \code{n}x\code{q} matrix of responses assumed to be drawn from a \code{q}-dimensional multivariate normal
  distribution.  \code{n} pertains to the sample size.}
  \item{x}{An \code{n}x\code{m} matrix of predictors (i.e., the design matrix). \code{m} pertains to the number of predictors. Do NOT include a column of
  1's if assuming an intercept.  This is controlled by the \code{int} argument.}
  \item{new.x}{A matrix of new values for which to approximate k-factors. This must be a matrix with \code{m} columns.}
  \item{int}{If \code{int = TRUE}, then an intercept is assumed.}
  \item{alpha}{The level chosen such that \code{1-alpha} is the confidence level.}
  \item{P}{The proportion of the population to be covered by this tolerance region.}
  \item{B}{The number of iterations used for the Monte Carlo algorithm which determines the tolerance factor.  The number of
  iterations should be at least as large as the default value of 1000.}
} \value{
  \code{mvregtol.region} returns a matrix where the first column is the k-factor, the next \code{q} columns are the estimated responses
  from the least squares fit, and the final \code{m} columns are the predictor values.  The first \code{n} rows of the matrix pertain to the raw data
  as specified by \code{y} and \code{x}.  If values for \code{new.x} are specified, then there is one additional row appended to this output for each
  row in the matrix \code{new.x}.
}

\references{
  Anderson, T. W. (2003) \emph{An Introduction to Multivariate Statistical Analysis}, Third Edition, Wiley.
  
  Krishnamoorthy, K. and Mathew, T. (2009), \emph{Statistical Tolerance Regions: Theory, Applications, and Computation}, Wiley.

  Krishnamoorthy, K. and Mondal, S. (2008), Tolerance Factors in Multiple and Multivariate Linear Regressions, 
  \emph{Communications in Statistics - Simulation and Computation}, \bold{37}, 546--559.
} 

\details{
  A basic sketch of how the algorithm works is as follows:
  
  (1) Generate independent chi-square random variables and Wishart random matrices.
  
  (2) Compute the eigenvalues of the randomly generated Wishart matrices.
  
  (3) Iterate the above steps to generate a set of \code{B} sample values such that the \code{100(1-alpha)}-th percentile is an approximate tolerance factor.
}


\examples{ 

## 95\%/95\% multivariate regression tolerance factors using
## a fertilizer data set presented in Anderson (2003, p. 374). 

grain <- c(40, 17, 9, 15, 6, 12, 5, 9)
straw <- c(53, 19, 10, 29, 13, 27, 19, 30)
fert <- c(24, 11, 5, 12, 7, 14, 11, 18)
Y <- cbind(grain, straw)
X <- cbind(fert)
new.x <- c(10, 15, 20)

set.seed(100)
out <- mvregtol.region(Y, X, new.x, int = TRUE, alpha = 0.05,
                       P = 0.95, B = 5000)
out
}

\keyword{file}




