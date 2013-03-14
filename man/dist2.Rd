\name{dist2}
\alias{dist2}

\title{
  Calculate an n-by-n matrix by applying a function to
  all pairs of columns of an m-by-n matrix.
}
\description{
  Calculate an n-by-n matrix by applying a function to
  all pairs of columns of an m-by-n matrix.
}
\usage{
  dist2(x, fun, diagonal=0)
}

\arguments{
  \item{x}{A matrix.}
  \item{fun}{A symmetric function of two arguments that may be columns of \code{x}.}
  \item{diagonal}{The value to be used for the diagonal elements of the
    resulting matrix.}
}
\details{
  With the default value of \code{fun}, this function calculates
  for each pair of columns of \code{x} the mean of the absolute values
  of their differences (which is proportional to the L1-norm of their
  difference). This is a distance metric.

  The implementation assumes that
  \code{fun(x[,i], x[,j])} can be evaluated for all pairs of \code{i}
  and \code{j} (see examples), and that
  \code{fun} is symmetric, i.e.
  \code{fun(a, b) = fun(b, a)}.
  \code{fun(a, a)} is not actually evaluated, instead the value of \code{diagonal}
  is used to fill the diagonal elements of the returned matrix.

  Note that \code{\link[stats:dist]{dist}} computes distances between rows of
  \code{x}, while this function computes relations between columns of
  \code{x} (see examples).
}
\value{
  A symmetric matrix of size \code{n x n}.
}
\author{
  Wolfgang Huber, James Reid
}
\examples{

  # example matrix
  z = matrix(1:15693, ncol=3)
  matL1 = dist2(z)
  matL2 = dist2(z, fun=function(a,b) sqrt(sum((a-b)^2, na.rm=TRUE)))

  euc = as.matrix(dist(t(z)))

  stopifnot(identical(dim(matL2), dim(euc)),
            all(euc==matL2))

}
\keyword{manip}








