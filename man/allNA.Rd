\name{allNA}
\alias{allNA}
\alias{anyNA}
\title{ A filter function to determine if all elements of a vector are NA. }
\description{
  \code{allNA} evaluates to \code{FALSE} if all elements of its
  argument are \code{NA}. \code{anyNA} evaluates to \code{FALSE} if any
  of the elements of its argument are \code{NA}.
}
\usage{
allNA(x)
anyNA(x)
}

\arguments{
  \item{x}{The vector to test.}
}

\value{
  \code{FALSE} if all elements of \code{x} are \code{NA}.
}

\author{R. Gentleman}

\seealso{\code{\link{poverA}}}

\examples{
  allNA(NA)
  allNA(1)
  anyNA(1)
  anyNA(NA)
}
\keyword{manip}
