\name{gapFilter}
\alias{gapFilter}
\title{ A filter to select genes based on there being a gap. }
\description{
  The \code{gapFilter} looks for genes that might usefully discriminate
  between two groups (possibly unknown at the time of filtering).
  To do this we look for a gap in the ordered expression values. The gap
  must come in the central portion (we exclude jumps in the initial
  \code{Prop} values or the final \code{Prop} values).
  Alternatively, if the IQR for the gene is large that will also pass
  our test and the gene will be selected.
}
\usage{
gapFilter(Gap, IQR, Prop, na.rm=TRUE, neg.rm=TRUE)
}
\arguments{
  \item{Gap}{The size of the gap required to pass the test. }
  \item{IQR}{The size of the IQR required to pass the test. }
  \item{Prop}{The proportion (or number) of samples to exclude at either
  end.}
  \item{na.rm}{If \code{TRUE} then \code{NA}'s will be removed before
    processing. }
  \item{na.rm}{If \code{TRUE} then elements of \code{x} that are negative
    will be removed before processing. }
}
\details{
  As stated above we are interested in 
}
\value{
 A function that returns either \code{TRUE} or \code{FALSE} depending on
 whether the vector supplied has a gap larger than \code{Gap} or an IQR
 (inter quartile range) larger than \code{IQR}. For computing the gap we
 want to exclude a proportion, \code{Prop} from either end of the sorted
 values. The reason for this requirement is that genes which differ in
 expression levels only for a few samples are not likely to be interesting.
}

\author{R. Gentleman }

\seealso{\code{\link{ttest}}, \code{\link{genefilter}} }

\examples{
 set.seed(256)
 x <- c(rnorm(10,100,3), rnorm(10, 100, 10))
 y <- x + c(rep(0,10), rep(100,10))
 tmp <- rbind(x,y) 
 Gfilter <- gapFilter(200, 100, 5)
 ffun <- filterfun(Gfilter)
 genefilter(tmp, ffun)
}
\keyword{manip}
