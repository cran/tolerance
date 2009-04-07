\name{acc.samp}
\title{Acceptance Sampling}
\alias{acc.samp}
\usage{
acc.samp(n, N, alpha = 0.05, P = 0.99, AQL)
}
\description{
  Provides an upper bound on the number of acceptable rejects or nonconformities in a process.  This is similar
  to a 1-sided upper tolerance bound for a hypergeometric random variable.
}

\arguments{
  \item{n}{The sample size to be drawn from the inventory.}
  \item{N}{The total inventory (or lot) size.}
  \item{alpha}{\code{1-alpha} is the confidence level for bounding the probability of accepting the inventory.}
  \item{P}{The proportion of items in the inventory which are to be accountable.}
  \item{AQL}{The acceptable quality level, which is the largest proportion of defects in a process considered
  acceptable.}

} \value{
  \code{acc.samp} returns a data frame with items:
  \item{acceptance.limit}{The number of items in the sample which may be unaccountable, yet still be able to 
  attain the desired confidence level \code{1-alpha}.}
  \item{lot.size}{The total inventory (or lot) size \code{N}.}
  \item{RQL}{The rejectable quality level.  This is the proportion of individual items in a sample one is willing
  to tolerate missing (i.e., this is \code{1-P}).}
  \item{confidence}{The confidence level \code{1-alpha}.}
  \item{AQL}{The acceptable quality level.  If the sampling were to be repeated numerous times as a process, then
  this quantity specifies the proportion of missing items considered acceptable from the process as a whole.}
  \item{sample.size}{The sample size drawn as specified by \code{n}.}
  \item{prod.risk}{The producer's risk.  This is the probability of rejecting an audit of a good inventory (also
  called the Type I error).  A good inventory can be rejected if an unfortunate random sample is selected (e.g.,
  most of the missing items happened to be selected for the audit).}
  \item{cons.risk}{The consumer's risk.  This is the probability of accepting an audit of a bad inventory (also
  called the Type II error).  A bad inventory can be accepted if a fortunate random sample (e.g., most of the missing
  items happend to not be selected for the audit).  \code{1-cons.risk} gives the actual confidence level of this 
  sampling plan.  If it is lower than the confidence level desired (i.e., because too small a sample size was
  specified), then a warning message will be displayed.}
} \seealso{
  \code{\link{Hypergeometric}}
}

\references{
  Montgomery, D. C. (2005), \emph{Introduction to Statistical Quality Control}, Fifth Edition, John Wiley \& Sons, Inc.
} 

\examples{ 

## A 90\%/90\% acceptance sampling plan for a sample of 450 
## drawn from a lot size of 960.

acc.samp(n = 450, N = 960, alpha = 0.10, P = 0.90, AQL = 0.01)
 }

\keyword{file}
