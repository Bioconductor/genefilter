\name{genescale}
\alias{genescale}
\alias{scaleVector}
\title{Scales a matrix or vector.}
\description{
  genescale() returns a scaled version of the input matrix m by applying
  the following formula to each column of the matrix:
     {    y[i] = ( x[i] - min(x) ) / ( max(x) - min(x) ) }

  scaleVector() returns a scaled version of the input vector v by applying
  the above formula to each element of the vector.
}
\usage{
genescale(m)

scaleVector(v)
}

\arguments{
  \item{m}{Input must be a matrix with numeric elements. }
  \item{v}{Input must be a vector with numeric elements. }

}
\details{
}
\value{
  Both functions return an identically sized, scaled version of the input.
}
\references{ }
\author{ M. Kajen }
\note{ These are helper functions used by genefinder(). }

\seealso{ \code{\link{genefinder}},\code{\link{scale}} }

\examples{
  v<-1:10
  scaleVector(v)

  m <- matrix(1:12, 4, 3)
  genescale(m)
}
\keyword{ manip }

