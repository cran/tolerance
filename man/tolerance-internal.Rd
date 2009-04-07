\name{tolerance-internal}
\alias{rwishart}

\title{Generation of a Random Wishart Matrix}
\description{
  Internal function for the package \code{tolerance}.
}
\usage{
rwishart(df, p)
}

\arguments{
  \item{df}{The degrees of freedom for the Wishart matrix to be generated.}
  \item{p}{The dimension of the random Wishart matrix.}
}

\details{
  This function is not intended to be called by the user.
}

\seealso{
  \code{\link{mvtol.region}}
}

\keyword{internal}
