\name{genefilter}
\alias{genefilter}

\title{A function to filter genes. }
\description{
   \code{genefilter} filters genes in the array \code{expr} using the
   filter functions in \code{flist}. It returns an array of logical
   values (suitable for subscripting) of the same length as there are
   rows in \code{expr}. For each row of \code{expr} the returned value
   is \code{TRUE} if the row passed all the filter functions. Otherwise
   it is set to \code{FALSE}.
}
\usage{
genefilter(expr, flist)
}
\arguments{
  \item{expr}{A matrix or dataframe that the filter functions will be
    applied to. }
  \item{flist}{A list of filter functions to apply to the array. }
}

\value{
  A logical vector of length equal to the number of rows of \code{expr}.
  The values in that vector indicate whether the corresponding row of
  \code{expr} passed the set of filter functions.
}
\author{R. Gentleman }
\seealso{\code{\link{genefilter}}, \code{\link{kOverA}} }

\examples{
   set.seed(-1)
   f1 <- kOverA(5, 10)
   flist <- filterfun(f1, allNA)
   exprA <- matrix(rnorm(1000, 10), nc=10)
   ans <- genefilter(exprA, flist)
}
\keyword{manip}

