\name{genescale}
\alias{genescale}
\title{Scales a matrix or vector.}
\description{
  \code{genescale} returns a scaled version of the input matrix m by applying
  the following formula to each column of the matrix:
     {    y[i] = ( x[i] - min(x) ) / ( max(x) - min(x) ) }
}
\usage{
genescale(m, axis=2, method=c("Z", "R"), na.rm=TRUE)
}

\arguments{
  \item{m}{Input a matrix or a vector with numeric elements. }
  \item{axis}{An integer indicating which axis of \code{m} to scale.}
  \item{method}{Either "Z" or "R", indicating whether a Z scaling or a
    range scaling should be performed.}
  \item{na.rm}{A boolean indicating whether \code{NA}'s should be
    removed.}
}
\details{
  Either the rows or columns of \code{m} are scaled. This is done either
  by subtracting the mean and dividing by the standard deviation ("Z")
  or by subtracing the minimum and dividing by the range.
}
\value{
  A scaled version of the input.
  If \code{m} is a \code{matrix} or a \code{dataframe} then the
  dimensions of the returned value agree with that of \code{m},
  in both cases the returned value is a \code{matrix}.
}

\author{ R. Gentleman }

\seealso{ \code{\link{genefinder}},\code{\link{scale}} }

\examples{
  m <- matrix(1:12, 4, 3)
  genescale(m)
}
\keyword{ manip }

