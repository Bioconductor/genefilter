\name{genescale}
\alias{genescale}
\alias{scaleVector}
\title{Scales a matrix or vector.}
\description{
  \code{genescale} returns a scaled version of the input matrix m by applying
  the following formula to each column of the matrix:
     {    y[i] = ( x[i] - min(x) ) / ( max(x) - min(x) ) }

  \code{scaleVector} returns a scaled version of the input vector v by applying
  the above formula to each element of the vector.
}
\usage{
genescale(m, axis)

scaleVector(v)
}

\arguments{
  \item{m}{Input a matrix or a vector with numeric elements. }
  \item{v}{Input must be a vector with numeric elements. }
  \item{axis}{An integer indicating which axis of \code{m} to apply
    \code{scaleVector} to. }
}
\details{
}
\value{
  Both functions return a scaled version of the input.
  If \code{m} is a \code{matrix} or a \code{dataframe} then the
  dimensions of the returned value agree with that of \code{m},
  in either case the returned value is a \code{matrix}.
}

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

