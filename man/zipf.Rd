\name{Zipf}
\alias{Zipf}
\alias{dzipf}
\alias{pzipf}
\alias{qzipf}
\alias{rzipf}
\title{
The Zipf Distribution
}
\description{
Density (mass), distribution function, quantile function, and random generation for the Zipf distribution with \code{N} categories and shape parameter \code{s}. 
}
\usage{
dzipf(x, s, N, log = FALSE)
pzipf(q, s, N, lower.tail = TRUE, log.p = FALSE)
qzipf(p, s, N, lower.tail = TRUE, log.p = FALSE)
rzipf(n, s, N)
}

\arguments{
    \item{x, q}{Vector of quantiles.}
    \item{p}{Vector of probabilities.}
    \item{n}{The number of observations.  If \code{length>1}, then the length is taken to be the number
    required.}
    \item{s}{The shape parameter, which must be greater than 0.}
    \item{N}{The number of categories, which must be integer-valued.}
    \item{log, log.p}{Logical vectors.  If \code{TRUE}, then the probabilities are given as \code{log(p)}.}
    \item{lower.tail}{Logical vector.  If \code{TRUE}, then probabilities are \eqn{P[X\le x]}, else \eqn{P[X>x]}.}
}
\details{
The Zipf distribution has mass
\deqn{p(x) = \frac{x^{-\lambda}}{\sum_{i=1}^{N}i^{-\lambda}},}{%
      p(x) = x^-\lambda/\sum_{i=1}^{N}i^(-\lambda),}	
where \eqn{x=1,\ldots,N}, \eqn{\lambda>0} is the shape parameter, and \code{N} is the number of distinct categories. Note that the Zipf distribution is just a special case of the Zipf-Mandelbrot distribution where the second shape parameter \code{b=0}.
\
}
\references{
Zipf, G. K. (1949), \emph{Human Behavior and the Principle of Least Effort}, Hafner.\\

Z\"{o}rnig, P. and Altmann, G. (1995), Unified Representation of Zipf Distributions, \emph{Computational Statistics and Data Analysis}, \bold{19}, 461--473. 
}
\value{
\code{dzipf} gives the density (mass), \code{pzipf} gives the distribution function, \code{qzipf} gives the quantile function, and \code{rzipf} generates random deviates. 
}
\seealso{
\code{\link{runif}} and \code{\link{.Random.seed}} about random number generation.
}
\examples{
## Randomly generated data from the Zipf distribution.

set.seed(100)
x <- rzipf(n = 500, s = 2, N = 100)
hist(x, main = "Randomly Generated Data", prob = TRUE)

x.1 <- sort(x)
y <- dzipf(x = x.1, s = 2, N = 100)
lines(x.1, y, col = 2, lwd = 2)

plot(x.1, pzipf(q = x.1, s = 2, N = 100), type = "l", 
     xlab = "x", ylab = "Cumulative Probabilities")

qzipf(p = 0.20, s = 2, N = 100, lower.tail = FALSE)
qzipf(p = 0.80, s = 2, N = 100)
}
\keyword{file}
