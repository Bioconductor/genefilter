\name{genefinder}
\alias{genefinder}
\title{ Finds genes that have similar patterns of expression. }
\description{
  Given a matrix of gene expressions, and the indexes of the genes of
  interest (in the matrix), genefinder() returns a list of distances
  between each gene of interest, and the entire matrix.
}
\usage{
genefinder(X, ilist, scale="none", method="euclidean")
}
\arguments{
  \item{X}{ A numeric matrix where columns represent patients, and rows
    represent genes. }
  \item{ilist}{ Vector of genes of interest. Contains indices of genes
    in matrix X. }
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
  specified in ilist. Each list entry is a vector containing the results
  of the distance comparison, with the gene specified in ilist, and
  every other row.
}
\author{ M. Kajen }

\seealso{\code{\link{genescale}} }

\examples{
m1 <- matrix (1:12, 4, 3)
v1 <- 1
genefinder (m1, v1)


v2 <- c(1,3)
genefinder (m1, v2)

genefinder (m1, v2, scale='range')

genefinder (m1, v2, method='manhattan')

m2 <- matrix (rnorm(100), 10, 10)
v3 <- c(2, 5, 6, 8)
genefinder (m2, v3, scale='zscore')
}
\keyword{ manip }

