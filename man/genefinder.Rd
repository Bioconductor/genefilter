\name{genefinder}
\alias{genefinder}
\title{ Finds genes that have similar patterns of expression. }
\description{
  Given a matrix of gene expressions, and the indices of the genes of
  interest (in the matrix), \code{genefinder} returns a list of the
  \code{numbResults} closest genes.
  The user can specify one of the standard distance measures listed
  below.
  The number of values to return can be specified. The return value is a
  list with two components:
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
  \item{scale}{ one of 'none', 'range', or 'zscore'. Scaling
    is carried out separately on each row.}
  \item{method}{ one of "euclidean", "maximum", "manhattan", "canberra", 
        "correlation", "binary".  }
}
\details{
  If the scale option is "range", then the input matrix is scaled using
  genescale(). If it is "zscore", then the input matrix is scaled using
  the 'scale' builtin with no arguments.

  The method option specifies the metric used for gene comparisons. The
  metric is applied, row by row, for each gene specified in \code{ilist}.

  The "correlation" option for the distance method will return a value
  equal to 1-correlation(x).
  
  See \code{\link{dist}} for a more detailed description of the distances.
}
\value{
  The returned value is a list containing an entry for each gene
  specified in ilist.  Each list entry contains an array of distances
  for that gene of interest.
}
\author{ J. Gentry and M. Kajen }

\seealso{\code{\link{genescale}} }

\examples{
set.seed(12345)
#create some fake expression profiles
m1 <- matrix (1:12, 4, 3)
v1 <- 1
nr <- 2
#find the 2 rows of m1 that are closest to row 1
genefinder (m1, v1, nr, method="euc")


v2 <- c(1,3)
genefinder (m1, v2, nr)

genefinder (m1, v2, nr, scale="range")

genefinder (m1, v2, nr, method="manhattan")

m2 <- matrix (rnorm(100), 10, 10)
v3 <- c(2, 5, 6, 8)
nr2 <- 6
genefinder (m2, v3, nr2, scale="zscore")

\testonly{
	m1 <- matrix(rnorm(1000),100,10)
	v1 <- c(3,5,8,42)
	nr2 <- 35
	genefinder(m1,v1,nr2,method="euclidean")
	genefinder(m1,v1,nr2,method="maximum")
	genefinder(m1,v1,nr2,method="canberra")
	genefinder(m1,v1,nr2,method="binary")
	genefinder(m1,v1,nr2,method="correlation")
	
	m2 <- matrix(rnorm(10000),1000,10)
	v1 <- c(1,100,563,872,921,3,52,95,235,333)
	nr <- 100
	genefinder(m2,v1,nr2,scale="zscore",method="euclidean")
	genefinder(m2,v1,nr2,scale="range", method="maximum")
	genefinder(m2,v1,nr2,scale="zscore", method="canberra")
	genefinder(m2,v1,nr2,scale="range", method="binary")
	genefinder(m2,v1,nr2,scale="zscore",method="correlation")
	}

}
\keyword{ manip }

