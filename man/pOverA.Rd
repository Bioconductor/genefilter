\name{pOverA}
\alias{pOverA}

\title{A filter function to filter according to the proportion of
  elements larger than A. }
\description{
  A function that returns a function with values for \code{A}, \code{p}
  and \code{na.rm} bound to the specified values. The function takes a
  single vector, \code{x}, as an argument.
  When the returned function is evaluated it returns \code{TRUE} if the
  proportion of values in \code{x} that are larger than \code{A} exceeds
  \code{p}. 
}
\usage{
pOverA(A=100, p=0.05, na.rm=TRUE)
}
\arguments{
  \item{A}{The value to be exceeded. }
  \item{p}{The proportion that need to exceed \code{A} for \code{TRUE}
    to be returned. }
  \item{na.rm}{ If \code{TRUE} then \code{NA}'s are removed. }
}

\value{
 \code{pOverA} returns a function with bindings for \code{A}, \code{p}
 and \code{na.rm}. This function evaluates to \code{TRUE} if the
  proportion of values in \code{x} that are larger than \code{A} exceeds
  \code{p}. 
}
\author{R. Gentleman}

\seealso{  \code{\link{sdom}} }

\examples{
  ff<- pOverA(10, p=.1)
  ff(1:20)
  ff(1:5)
}
\keyword{manip }
