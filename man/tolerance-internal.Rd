\name{tolerance-internal}
\alias{rwishart}
\alias{zeta.fun}

\title{Internal Functions}
\description{
Internal functions for the package \code{tolerance}.
}
\usage{
rwishart(df, p)
zeta.fun(x)
}

\arguments{
  \item{df}{The degrees of freedom for the Wishart matrix to be generated.}
  \item{p}{The dimension of the random Wishart matrix.}
  \item{x}{A vector or matrix whose real values must be greater than or equal to 1.}
}

\details{
These functions are not intended to be called by the user.  \code{rwishart} performs random generation of Wishart matrices.  \code{zeta.fun} is a condensed version of the Riemann's zeta function given in the VGAM package.  Please use that reference if looking to directly implement Riemann's zeta function.  The function we have included is done so out of convenience.  
}


\references{
  Yee, T. (2010), The VGAM Package for Categorical Data Analysis, \emph{Journal of Statistical Software}, \bold{32}, 1--34.
} 


\seealso{
  \code{\link{mvtol.region}}
}

\keyword{internal}
