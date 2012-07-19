\name{tolerance-internal}
\alias{rFUN}
\alias{rwishart}
\alias{zeta.fun}

\title{Internal Functions}
\description{
Internal functions for the package \code{tolerance}.
}
\usage{
rFUN(FUN, r1 = "1", r2 = "2")
rwishart(df, p)
zeta.fun(x)
}

\arguments{
  \item{FUN}{A function that takes exactly two arguments: \code{r1} and \code{r2}.}
  \item{r1, r2}{The two arguments to \code{FUN}.}
  \item{df}{The degrees of freedom for the Wishart matrix to be generated.}
  \item{p}{The dimension of the random Wishart matrix.}
  \item{x}{A vector or matrix whose real values must be greater than or equal to 1.}
}

\details{
These functions are not intended to be called by the user.  \code{rFUN} is a way to extract the character expression for \code{FUN} with the two arguments \code{r1} and \code{r2}. \code{rwishart} performs random generation of Wishart matrices.  \code{zeta.fun} is a condensed version of the Riemann's zeta function given in the VGAM package.  Please use that reference if looking to directly implement Riemann's zeta function.  The function we have included is done so out of convenience.  
}


\references{
  Yee, T. (2010), The VGAM Package for Categorical Data Analysis, \emph{Journal of Statistical Software}, \bold{32}, 1--34.
} 


\seealso{
  \code{\link{mvtol.region}}
}

\keyword{internal}
