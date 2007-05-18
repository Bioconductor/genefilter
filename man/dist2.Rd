\name{dist2}
\alias{dist2}

\title{
  Calculate an n-by-n matrix by applying a function to 
  pairs of columns of an m-by-n matrix.
}
\description{
  Calculate an n-by-n matrix by applying a function to 
  pairs of columns of an m-by-n matrix.
}
\usage{
  dist2(x, fun=function(a,b) mad(a-b))
}

\arguments{
  \item{x}{A matrix, or indeed any object for which \code{ncol(x)}
    and \code{x[,j]} return useful results.}
  \item{fun}{A symmetric function of two arguments that may be
    columns of \code{x}.}
}
\details{
  With the default value of \code{fun}, this calculates a matrix of the
  \code{\link{mad}} of all pairwise differences of columns of \code{x}.
  This may be considered a measure of distance.

  The implementation assumes that \code{fun} is symmetric,
  \code{fun(a,b)=fun(b,a)}.
  
  A use for this function is the detection of outlier arrays in a
  microarray experiment. Assume that each column of \code{x} can be
  decomposed as
  $z+\beta+\epsilon$, where $z$ is a fixed vector
  (the same for all columns), $\epsilon$ is vector of
  \code{nrow{x}} i.i.d. random numbers, and $\beta$ is an arbitrary
  vector whose majority of entries are negligibly small (i.e. close to
  zero). In other words, $z$ the probe effects,  $\epsilon$ measurement
  noise and $\beta$ differential expression effects. Under this
  assumption, all entries of the resulting distance matrix should be the
  same, namely $\sqrt(2)$ times the standard deviation of $\epsilon$.
  Arrays whose  distance matrix entries are way different give cause
  for suspicion.
}  
\value{
  A symmetric matrix of size \code{n x n}.
  Diagonal elements are \code{NA}.
}

\examples{
  z = matrix(rnorm(15693), ncol=3)
  dist2(z)
}
\keyword{manip}








