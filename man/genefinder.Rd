\name{genefinder}
\alias{genefinder}
\title{ Finds genes that have similar patterns of expression. }
\description{
  Given a matrix of gene expressions, and the indices of the genes of
  interest (in the matrix), genefinder() returns a list of the X closest
  genes (measured through the desired distance method) to the genes
  of interest (where X is the number of desired results returned) and
  their distances.
}
\usage{
genefinder(X, ilist, numResults=25, scale="none", method="euclidean")
}
\arguments{
  \item{X}{ A numeric matrix where columns represent patients, and rows
    represent genes. }
  \item{ilist}{ Vector of genes of interest. Contains indices of genes
    in matrix X. }
  \item{numResults}{Number of results to display, starting from the least
    distance to the greatest. }
  \item{scale}{ one of 'none', 'range', or 'zscore'. Selecting 'range'
    causes the input matrix to be scaled before other computations }
  \item{method}{ one of "euclidean", "maximum", "manhattan", "canberra", 
        "binary".  }
}
\details{
  If the scale option is 'range', then the input matrix is scaled using
  genescale(). If it is 'zscore', then the input matrix is scaled using
  the 'scale' builtin with no arguments.

  The method option specifies the metric used for gene comparisons. The
  metric is applied, column by column, for each gene specified in the ilist.

}
\value{
  The returned value is a list containing an entry for each gene
  specified in ilist.  Each list entry contains an array of distances
  for that gene of interest.
}
\author{ M. Kajen }

\seealso{\code{\link{genescale}} }

\examples{
m1 <- matrix (1:12, 4, 3)
v1 <- 1
nr <- 2
genefinder (m1, v1, nr)


v2 <- c(1,3)
genefinder (m1, v2, nr)

genefinder (m1, v2, nr, scale='range')

genefinder (m1, v2, nr, method='manhattan')

m2 <- matrix (rnorm(100), 10, 10)
v3 <- c(2, 5, 6, 8)
nr2 <- 6
genefinder (m2, v3, nr2, scale='zscore')
}
\keyword{ manip }

