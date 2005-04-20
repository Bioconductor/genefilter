\name{rowSds}
\alias{rowSds}
\title{Row standard deviation of a numeric array}
\description{
  Row standard deviation of a numeric array
}
\usage{
rowSds(x, ...)
}
\arguments{
  \item{x}{An array of two or more dimensions, containing numeric,
          complex, integer or logical values, or a numeric data frame.}
  \item{...}{Further arguments that get passed on to
    \code{\link{rowMeans}} and \code{\link{rowSums}}.}
}

\value{
  A numeric or complex array of suitable size, or a vector if the
  result is one-dimensional.  The `dimnames' (or `names' for a
  vector result) are taken from the original array.
}
\details{This a convenience function, the main work is done in
  \code{\link{rowMeans}} and \code{\link{rowSums}}. See the function
  definition, it is very simple.
}
\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}
\seealso{\code{\link{rowMeans}} and \code{\link{rowSums}}}

\examples{
   a = matrix(rnorm(1e4), nrow=10)
   rowSds(a)
}

\keyword{array}
\keyword{manip}
