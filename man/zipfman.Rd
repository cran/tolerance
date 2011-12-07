\name{ZipfMandelbrot}
\alias{ZipfMandelbrot}
\alias{dzipfman}
\alias{pzipfman}
\alias{qzipfman}
\alias{rzipfman}
\title{
The Zipf-Mandelbrot Distribution
}
\description{
Density (mass), distribution function, quantile function, and random generation for the Zipf-Mandelbrot distribution with \code{N} categories and shape parameters \code{s} and \code{b}. 
}
\usage{dzipfman(x, s, b, N, log = FALSE)
pzipfman(q, s, b, N, lower.tail = TRUE, log.p = FALSE)
qzipfman(p, s, b, N, lower.tail = TRUE, log.p = FALSE)
rzipfman(n, s, b, N)
}

\arguments{
    \item{x, q}{Vector of quantiles.}
    \item{p}{Vector of probabilities.}
    \item{n}{The number of observations.  If \code{length>1}, then the length is taken to be the number
    required.}
    \item{s, b}{The shape parameters, both of which must be greater than 0.}
    \item{N}{The number of categories, which must be integer-valued.}
    \item{log, log.p}{Logical vectors.  If \code{TRUE}, then the probabilities are given as \code{log(p)}.}
    \item{lower.tail}{Logical vector.  If \code{TRUE}, then probabilities are \eqn{P[X\le x]}, else \eqn{P[X>x]}.}
}
\details{
The Zipf-Mandelbrot distribution has mass
\deqn{p(x) = \frac{(x + \delta)^{-\lambda}}{\sum_{i=1}^{N}(i + \delta)^{-\lambda}},}{%
      p(x) = (x + \delta)^-\lambda/\sum_{i=1}^{N}(i + \delta)^(-\lambda),}	
where \eqn{x=1,\ldots,N}, \eqn{\lambda,\delta>0} are shape parameters, and \code{N} is the number of distinct categories. Note that the Zipf-Mandelbrot distribution is just a generalization of the Zipf distribution where the second shape parameter \code{b=0}.
}
\references{
Mandelbrot, B. B. (1965), Information Theory and Psycholinguistics. In B. B. Wolman and E. Nagel, editors. \emph{Scientific Psychology}, Basic Books.\\

Zipf, G. K. (1949), \emph{Human Behavior and the Principle of Least Effort}, Hafner.\\

Z\"{o}rnig, P. and Altmann, G. (1995), Unified Representation of Zipf Distributions, \emph{Computational Statistics and Data Analysis}, \bold{19}, 461--473. 
}
\value{
\code{dzipfman} gives the density (mass), \code{pzipfman} gives the distribution function, \code{qzipfman} gives the quantile function, and \code{rzipfman} generates random deviates. 
}
\seealso{
\code{\link{runif}} and \code{\link{.Random.seed}} about random number generation.
}
\examples{
## Randomly generated data from the Zipf-Mandelbrot distribution.

set.seed(100)
x <- rzipfman(n = 500, s = 2, b = 3, N = 100)
hist(x, main = "Randomly Generated Data", prob = TRUE)

x.1 <- sort(x)
y <- dzipfman(x = x.1, s = 2, b = 3, N = 100)
lines(x.1, y, col = 2, lwd = 2)

plot(x.1, pzipfman(q = x.1, s = 2, b = 3, N = 100), type = "l", 
     xlab = "x", ylab = "Cumulative Probabilities")

qzipfman(p = 0.20, s = 2, b = 3, N = 100, lower.tail = FALSE)
qzipfman(p = 0.80, s = 2, b = 3, N = 100)
}
\keyword{file}
