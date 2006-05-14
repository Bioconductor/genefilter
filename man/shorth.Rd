\name{shorth}
\alias{shorth}
\title{Midpoint of the shorth}
\description{Calculate the midpoint of the shorth}
\usage{shorth(x, na.rm=FALSE, tie.action="mean", tie.limit=0.05)}
\arguments{
  \item{x}{Numeric}
  \item{na.rm}{Logical. If \code{TRUE}, then \code{NA} values are ignored.}
  \item{tie.action}{Character scalar. See details.}
  \item{tie.limit}{Numeric scalar. See details.}
}

\details{The shorth is the shortest interval that covers half of the
  values in \code{x}. Its midpoint is the median of the data values
  that lie in the shorth and can be used as an
  estimator of the mode of the distribution that generated \code{x}.
  
  If there are multiple shortest intervals,
  the action specified in \code{ties.action} is applied.
  Allowed values are \code{mean} (the default), \code{max} and \code{min}.
  For \code{mean}, the average value is considered; however, an error is
  generated if the start indices of the different shortest intervals
  differ by more than the fraction \code{tie.limit} of \code{length(x)}.
  For \code{min} and \code{max}, the left-most or right-most, respectively, of
  the multiple shortest intervals is considered.
}

\value{The midpoint of the short of \code{x}.}
\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}, Ligia Pedroso Bras}
\seealso{\code{\link{half.range.mode}}}

\examples{
 
  x <- c(rnorm(500), runif(500) * 10)
  res <- c(mean=mean(x), median=median(x), shorth=shorth(x))
  if(interactive()) {
    colors = c(1,2,4)
    hist(x, 40, col="orange")
    abline(v=res, col=colors, lwd=3)
    legend(6, 100, names(res), col=colors, lwd=3) 
  }
}
\keyword{arith}

