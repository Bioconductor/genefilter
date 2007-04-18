\name{shorth}
\alias{shorth}
\title{A location estimator based on the shorth}
\description{A location estimator based on the shorth}
\usage{shorth(x, na.rm=FALSE, tie.action="mean", tie.limit=0.05)}
\arguments{
  \item{x}{Numeric}
  \item{na.rm}{Logical. If \code{TRUE}, then non-finite (according to
    \code{\link{is.finite}}) values in \code{x} are ignored. Otherwise,
    presence of non-finite or \code{NA} values will lead to an error message.}
  \item{tie.action}{Character scalar. See details.}
  \item{tie.limit}{Numeric scalar. See details.}
}

\details{The shorth is the shortest interval that covers half of the
  values in \code{x}. This function calculates the mean of the \code{x}
  values that lie in the shorth, and this can be used as a robust
  estimator of location.

  G. Sawitzki points out that, as an estimator of location, the quantity
  computed here has an asymptotic rate of only $n^{-1/3}$ and a
  complicated limiting distribution. 
  Hence, the function should be considered only for large samples, such
  as from microarray experiments, where the rate of convergence is not
  the prior concern.

  Ties: if there are multiple shortest intervals,
  the action specified in \code{ties.action} is applied.
  Allowed values are \code{mean} (the default), \code{max} and \code{min}.
  For \code{mean}, the average value is considered; however, an error is
  generated if the start indices of the different shortest intervals
  differ by more than the fraction \code{tie.limit} of \code{length(x)}.
  For \code{min} and \code{max}, the left-most or right-most, respectively, of
  the multiple shortest intervals is considered.

  See \code{\link{half.range.mode}} for an iterative version
  that refines the estimate iteratively and has a builtin bootstrapping option.
}

\value{The mean of the \code{x} values that lie in the shorth.}

\references{
  \itemize{
    \item G Sawitzki, \dQuote{Diagnostic Plots for One-Dimensional Data.}
    Presented at the 24th meeting of the Arbeitsgruppe "Computational
    Statistics", Internationale Biometrische Gesellschaft, Reisensburg (1992).
    
    \item DF Andrews, \dQuote{Robust Estimates of Location.}
    Princeton University Press (1972).

    \item R Gruebel, \dQuote{The Length of the Shorth.} Annals of
    Statistics 16, 2:619-628 (1988).
  }
}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}, Ligia Pedroso Bras}
\seealso{\code{\link{half.range.mode}}}

\examples{
 
  x = c(rnorm(500), runif(500) * 10)
  methods = c("mean", "median", "shorth", "half.range.mode")
  ests = sapply(methods, function(m) get(m)(x))

  if(interactive()) {
    colors = 1:4
    hist(x, 40, col="orange")
    abline(v=ests, col=colors, lwd=3, lty=1:2)
    legend(5, 100, names(ests), col=colors, lwd=3, lty=1:2) 
  }
}
\keyword{arith}

