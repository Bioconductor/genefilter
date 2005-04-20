\name{shorth}
\alias{shorth}
\title{Midpoint of the shorth}
\description{Calculate the midpoint of the shorth}
\usage{shorth(x, na.rm=FALSE)}
\arguments{
  \item{x}{Numeric}
  \item{na.rm}{Logical}
}

\details{The shorth is the shortest interval that covers half of the
  values in \code{x}. Its midpoint is the median of the data values
  that lie in the shorth and can be used as an
  estimator of the mode of the distribution that generated \code{x}.}

\value{The midpoint of the short of \code{x}.}
\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}

\examples{
  x <- c(rnorm(100), runif(100) * 10)
  mean(x)
  median(x)
  shorth(x)
}
\keyword{arith}

