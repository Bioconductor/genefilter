\name{ttest}
\alias{ttest}
\title{A filter function for a t.test }
\description{
  \code{ttest} returns a function of one argument with bindings for
  \code{cov} and \code{p}.
  The function, when evaluated, performs an ANOVA using \code{cov} as
  the covariate. It returns \code{TRUE} if the p value for a difference
  in means is less than \code{p}. 
}
\usage{
ttest(m, p=0.05, na.rm=TRUE)
}
\arguments{
  \item{m}{If \code{m} is of length one then it is assumed that elements
  one through \code{m} of \code{x} will be one group. Otherwise \code{m}
  is presumed to be the same length as \code{x} and constitutes the
  groups.}
  \item{p}{ The p-value for the test. }
  \item{na.rm}{If set to \code{TRUE} any \code{NA}'s will be removed. }
}

\value{
 \code{ttest} returns a function with bindings for \code{m} and
 \code{p} that will perform a \code{\link{t.test}}.
}
\author{R. Gentleman }
\seealso{\code{\link{kOverA}}, \code{\link{Anova}} }

\examples{
  dat <- c(rep(1,5),rep(2,5))
  set.seed(5)
  y <- rnorm(10)
  af <- ttest(dat, .01)
  af(y)
  af2 <- ttest(5, .01)
  af2(y)
  y[8] <- NA
  af(y)
  af2(y)
  y[1:5] <- y[1:5]+10
  af(y)
}
\keyword{manip}
