\name{nptol.int}
\title{Nonparametric Tolerance Intervals}
\alias{nptol.int}
\usage{
nptol.int(x, alpha = 0.05, P = 0.99, side = 1,
          method = c("WILKS", "WALD", "HM"), upper = NULL,
          lower = NULL)
}
\description{
  Provides 1-sided or 2-sided nonparametric (i.e., distribution-free) tolerance intervals for any continuous
  data set.
}

\arguments{
  \item{x}{A vector of data which no distributional assumptions are made.  The data is only assumed to come from a continuous
  distribution.}
  \item{alpha}{The level chosen such that \code{1-alpha} is the confidence level.}
  \item{P}{The proportion of the population to be covered by this tolerance interval.}
  \item{side}{Whether a 1-sided or 2-sided tolerance interval is required (determined by \code{side = 1} or \code{side = 2},
  respectively).}
  \item{method}{The method for determining which indices of the ordered observations will be used for the tolerance intervals.  
  \code{"WILKS"} is the Wilks method, which produces tolerance bounds symmetric about the observed center of the residuals by using
  the beta distribution.  \code{"WALD"} is the Wald method, which produces (possibly) multiple tolerance bounds for \code{side = 2} (each
  having at least the specified confidence level), but is the same as \code{method = "WILKS"} for \code{side = 1}.  \code{"HM"} is
  the Hahn-Meeker method, which is based on the binomial distribution, but the upper and lower bounds may exceed the minimum and maximum
  of the sample data.  For \code{side = 2}, this method will yield two intervals if an odd number of observations are to be trimmed from each side.}
  \item{upper}{The upper bound of the data.  When \code{NULL}, then the maximum of \code{x} is used.}
  \item{lower}{The lower bound of the data.  When \code{NULL}, then the minimum of \code{x} is used.}
} \value{
  \code{nptol.int} returns a data frame with items:
  \item{alpha}{The specified significance level.}
  \item{P}{The proportion of the population covered by this tolerance interval.}
  \item{1-sided.lower}{The 1-sided lower tolerance bound.  This is given only if \code{side = 1}.}
  \item{1-sided.upper}{The 1-sided upper tolerance bound.  This is given only if \code{side = 1}.}
  \item{2-sided.lower}{The 2-sided lower tolerance bound.  This is given only if \code{side = 2}.}
  \item{2-sided.upper}{The 2-sided upper tolerance bound.  This is given only if \code{side = 2}.}
} \seealso{
  \code{\link{distfree.est}}, \code{\link{npregtol.int}}
}

\references{
  Bury, K. (1999), \emph{Statistical Distributions in Engineering}, Cambridge University Press.

  Hahn, G. J. and Meeker, W. Q. (1991), \emph{Statistical Intervals: A Guide for Practitioners}, Wiley-Interscience.

  Wald, A. (1943), An Extension of Wilks' Method for Setting Tolerance Limits, \emph{The Annals of Mathematical Statistics},
  \bold{14}, 45--55.
  
  Wilks, S. S. (1941), Determination of Sample Sizes for Setting Tolerance Limits, \emph{The Annals of Mathematical Statistics},
  \bold{12}, 91--96.
} 

\examples{ 
## 90\%/95\% 2-sided nonparametric tolerance intervals for a 
## sample of size 20. 

set.seed(100)
x <- rlogis(20, 5, 1)
out <- nptol.int(x = x, alpha = 0.10, P = 0.95, side = 1,
                 method = "WILKS", upper = NULL, lower = NULL)
out

plottol(out, x, plot.type = "both", side = "two", x.lab = "X")
}

\keyword{file}




