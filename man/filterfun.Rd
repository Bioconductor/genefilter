\name{filterfun}
\alias{filterfun}
\title{Creates a first FALSE exiting function from the list of
  filter functions it is given. }
\description{
  This function creates a function that takes a single argument. The
  filtering functions are bound in the environment of the returned
  function and are applied sequentially to the argument of the returned
  function. When the first filter function evaluates to \code{FALSE} the
  function returns \code{FALSE} otherwise it returns \code{TRUE}.
}
\usage{
filterfun(...)
}
\arguments{
  \item{...}{Filtering functions. }
}

\value{
 \code{filterfun} returns a function that takes a single argument. It
 binds the filter functions given to it in the environment of the
 returned function. These functions are applied sequentially (in the
 order they were given to \code{filterfun}). The function returns
 \code{FALSE} (and exits) when the first filter function returns
 \code{FALSE} otherwise it returns \code{TRUE}.
}
\author{R. Gentleman }
\seealso{\code{\link{genefilter}} }
\examples{
 set.seed(333)
 x <- matrix(rnorm(100,2,1),nc=10)
 cvfun <- cv(.5,2.5)
 ffun <- filterfun(cvfun)
 which <- genefilter(x, ffun)
}
\keyword{manip}

