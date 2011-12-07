\name{Zeta}
\alias{Zeta}
\alias{dzeta}
\alias{pzeta}
\alias{qzeta}
\alias{rzeta}
\title{
The Zeta Distribution
}
\description{
Density (mass), distribution function, quantile function and random generation for the zeta distribution with parameter \code{s}.
}
\usage{
dzeta(x, s, log = FALSE)
pzeta(q, s, lower.tail = TRUE, log.p = FALSE)
qzeta(p, s, lower.tail = TRUE, log.p = FALSE)
rzeta(n, s)
}

\arguments{
    \item{x, q}{Vector of quantiles.}
    \item{p}{Vector of probabilities.}
    \item{n}{The number of observations.  If \code{length>1}, then the length is taken to be the number
    required.}
    \item{s}{The shape parameter, which must be strictly greater than 1.}
    \item{log, log.p}{Logical vectors.  If \code{TRUE}, then the probabilities are given as \code{log(p)}.}
    \item{lower.tail}{Logical vector.  If \code{TRUE}, then probabilities are \eqn{P[X\le x]}, else \eqn{P[X>x]}.}
}
\details{
The zeta distribution has mass
\deqn{p(x) = \frac{x^{-\lambda}}{\zeta(\lambda)},}{%
      p(x) = x^-\lambda/\zeta(\lambda),}	
where \eqn{x=1,2,\ldots}, \eqn{\lambda>1} is the shape parameter, and \eqn{\zeta()} is the Riemann zeta function given by:
\deqn{\zeta(t) = \sum_{i=1}^{\infty}\frac{1}{i^{t}}<\infty.}{%
			\zeta(t) = \sum_{i=1}^\infty 1/i^t<\infty.}
Note that the zeta distribution is just a special case of the Zipf distribution where \code{s>1} and \code{N} goes to infinity.
}

\references{
Zipf, G. K. (1949), \emph{Human Behavior and the Principle of Least Effort}, Hafner.\\

Z\"{o}rnig, P. and Altmann, G. (1995), Unified Representation of Zipf Distributions, \emph{Computational Statistics and Data Analysis}, \bold{19}, 461--473. 
}
\value{
\code{dzeta} gives the density (mass), \code{pzeta} gives the distribution function, \code{qzeta} gives the quantile function, and \code{rzeta} generates random deviates. 
}
\seealso{
\code{\link{runif}} and \code{\link{.Random.seed}} about random number generation.
}
\examples{
## Randomly generated data from the zeta distribution.

set.seed(100)
x <- rzeta(n = 100, s = 1.3)
hist(x, main = "Randomly Generated Data", prob = TRUE)

x.1 <- sort(x)
y <- dzeta(x = x.1, s = 1.3)
lines(x.1, y, col = 2, lwd = 2)

plot(x.1, pzeta(q = x.1, s = 1.3), type = "l", 
     xlab = "x", ylab = "Cumulative Probabilities")

qzeta(p = 0.20, s = 1.3, lower.tail = FALSE)
qzeta(p = 0.80, s = 1.3)
}
\keyword{file}
