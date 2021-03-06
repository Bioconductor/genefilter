\name{half.range.mode}
\alias{half.range.mode}
\title{Mode estimation for continuous data}

\description{
  For data assumed to be drawn from a unimodal, continuous distribution,
  the mode is estimated by the \dQuote{half-range} method. Bootstrap resampling
  for variance reduction may optionally be used.
}

\usage{
half.range.mode(data, B, B.sample, beta = 0.5, diag = FALSE)
}

\arguments{
  \item{data}{A numeric vector of data from which to estimate the mode.}
  \item{B}{
    Optionally, the number of bootstrap resampling rounds to use. Note
    that \code{B = 1} resamples 1 time, whereas omitting \code{B}
    uses \code{data} as is, without resampling.
  }
  \item{B.sample}{
    If bootstrap resampling is requested, the size of the bootstrap
    samples drawn from \code{data}. Default is to use a sample which is
    the same size as \code{data}. For large data sets, this may be slow
    and unnecessary.
  }
  \item{beta}{
    The fraction of the remaining range to use at each iteration.
  }
  \item{diag}{
    Print extensive diagnostics. For internal testing only... best left
    \code{FALSE}.
  }
}

\details{
  Briefly, the mode estimator is computed by iteratively identifying
  densest half ranges. (Other fractions of the current range can be
  requested by setting \code{beta} to something other than 0.5.) A densest half
  range is an interval whose width equals half the current range, and
  which contains the maximal number of observations. The subset of
  observations falling in the selected densest half range is then used to compute
  a new range, and the procedure is iterated. See the references for
  details.

  If bootstrapping is requested, \code{B} half-range mode estimates are
  computed for \code{B} bootstrap samples, and their average is returned
  as the final estimate.
}

\value{
  The mode estimate.
}

\references{
  \itemize{
    \item DR Bickel, \dQuote{Robust estimators of the mode and skewness of
    continuous data.} \emph{Computational Statistics & Data Analysis}
    39:153-163 (2002).

    \item SB Hedges and P Shah, \dQuote{Comparison of
    mode estimation methods and application in molecular clock analysis.} \emph{BMC
    Bioinformatics} 4:31-41 (2003).
  }
}

\author{Richard Bourgon <bourgon@stat.berkeley.edu>}
\seealso{\code{\link{shorth}}}
\keyword{univar}
\keyword{robust}

\examples{
## A single normal-mixture data set

x <- c( rnorm(10000), rnorm(2000, mean = 3) )
M <- half.range.mode( x )
M.bs <- half.range.mode( x, B = 100 )

if(interactive()){
hist( x, breaks = 40 )
abline( v = c( M, M.bs ), col = "red", lty = 1:2 )
legend(
       1.5, par("usr")[4],
       c( "Half-range mode", "With bootstrapping (B = 100)" ),
       lwd = 1, lty = 1:2, cex = .8, col = "red"
       )
}

# Sampling distribution, with and without bootstrapping

X <- rbind(
           matrix( rnorm(1000 * 100), ncol = 100 ),
           matrix( rnorm(200 * 100, mean = 3), ncol = 100 )
           )
M.list <- list(
               Simple = apply( X, 2, half.range.mode ),
               BS = apply( X, 2, half.range.mode, B = 100 )
               )

if(interactive()){
boxplot( M.list, main = "Effect of bootstrapping" )
abline( h = 0, col = "red" )
}
}
