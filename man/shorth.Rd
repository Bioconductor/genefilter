\name{shorth}
\alias{shorth}
\title{Midpoint of the shorth}
\description{Calculate the midpoint of the shorth}
\usage{shorth(x, na.rm=FALSE, tieLimit=0.05)}
\arguments{
  \item{x}{Numeric}
  \item{na.rm}{Logical}
  \item{tieLimit}{Numeric}
}

\details{The shorth is the shortest interval that covers half of the
  values in \code{x}. Its midpoint is the median of the data values
  that lie in the shorth and can be used as an
  estimator of the mode of the distribution that generated \code{x}.

  An error is generated if there are multiple shortest intervals whose
  start indices differ by more than a fraction \code{tieLimit} of the
  length of \code{x}.
}

\value{The midpoint of the short of \code{x}.}
\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}

\examples{
  x <- c(rnorm(100), runif(100) * 10)
  mean(x)
  median(x)
  shorth(x)
}
\keyword{arith}

