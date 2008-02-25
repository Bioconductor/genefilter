\name{genefilter}
\alias{genefilter}
\title{A function to filter genes.}
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
  \item{expr}{A \code{matrix} or \code{ExpressionSet} that the
    filter functions will be applied to.}
  \item{flist}{A \code{list} of filter functions to apply to the array.}
}
\details{
   This package uses a very simple but powerful protocol for
   \emph{filtering} genes. The user simply constructs any number of
   tests that they want to apply. A test is simply a function (as
   constructed using one of the many helper functions in this package)
   that returns \code{TRUE} if the gene of interest passes the test (or
   filter) and \code{FALSE} if the gene of interest fails.

   The benefit of this approach is that each test is constructed
   individually (and can be tested individually). The tests are then
   applied sequentially to each gene. The function returns a logical
   vector indicating whether the gene passed all tests functions or
   failed at least one of them.
}
\value{
  A logical \code{vector} of length equal to the number of rows of
  \code{expr}. The values in that \code{vector} indicate whether the
  corresponding row of \code{expr} passed the set of filter functions.
}
\author{R. Gentleman}
\seealso{\code{\link{genefilter}}, \code{\link{kOverA}}}
\examples{
   set.seed(-1)
   f1 <- kOverA(5, 10)
   flist <- filterfun(f1, allNA)
   exprA <- matrix(rnorm(1000, 10), ncol = 10)
   ans <- genefilter(exprA, flist)
}
\keyword{manip}
