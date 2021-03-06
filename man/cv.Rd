\name{cv}
\alias{cv}
\title{A filter function for the coefficient of variation.}
\description{
  \code{cv} returns a function with values for \code{a} and \code{b}
  bound. This function takes a single argument. It computes the
  coefficient of variation for the input vector and returns \code{TRUE} if
  the coefficient of variation is between \code{a} and
  \code{b}. Otherwise it returns \code{FALSE}
}
\usage{
cv(a=1, b=Inf, na.rm=TRUE)
}
\arguments{
  \item{a}{The lower bound for the cv. }
  \item{b}{The upper bound for the cv. }
  \item{na.rm}{If set to \code{TRUE} any \code{NA}'s will be removed. }
}
\details{
  The coefficient of variation is the standard deviation divided by the 
  absolute value of the mean.
}
\value{
 It returns a function of one argument. The function has an environment
 with bindings for \code{a} and \code{b}.
}
\author{R. Gentleman }
\seealso{\code{\link{pOverA}}, \code{\link{kOverA}} }

\examples{
  set.seed(-3)
  cvfun <- cv(1,10)
  cvfun(rnorm(10,10))
  cvfun(rnorm(10))
}
\keyword{manip}

