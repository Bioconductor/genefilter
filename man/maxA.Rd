\name{maxA}
\alias{maxA}
\title{ A filter function to filter according to the maximum. }
\description{
  \code{maxA} returns a function with the parameter \code{A} bound.
  The returned function evaluates to \code{TRUE} if any element of its
  argument is larger than \code{A}.
}
\usage{
maxA(A=75, na.rm=TRUE)
}
\arguments{
  \item{A}{The value that at least one element must exceed. }
  \item{na.rm}{If \code{TRUE} then \code{NA}'s are removed. }
}

\value{
 \code{maxA} returns a function with an environment containing a binding
 for \code{A}.
 
}

\author{R. Gentleman }

\seealso{\code{\link{pOverA}} }

\examples{
   ff <- maxA(30)
   ff(1:10)
   ff(28:31)
}
\keyword{manip}

