\name{kOverA}
\alias{kOverA}
\title{A filter function for k elements larger than A. }
\description{
  \code{kOverA} returns a filter function with bindings for \code{k} and
  \code{A}. This function evaluates to \code{TRUE} if at least \code{k}
  of the arguments elements are larger than \code{A}.
}
\usage{
kOverA(k, A=100, na.rm=TRUE)
}
\arguments{
  \item{A}{The value you want to exceed. }
  \item{k}{The number of elements that have to exceed A.}
  \item{na.rm}{If set to \code{TRUE} any \code{NA}'s will be removed. }
}

\value{
 A function with bindings for \code{A} and \code{k}.
}

\author{R. Gentleman}

\seealso{\code{\link{pOverA}}}

\examples{
   fg <- kOverA(5, 100)
   fg(90:100)
   fg(98:110)
}
\keyword{manip}
