\name{findLargest}
\alias{findLargest}
\title{Find the Entrez Gene ID corresponding to the largest statistic}
\description{
  Most microarrays have multiple probes per gene (Entrez). This function
  finds all replicates, and then selects the one with the largest value
  of the test statistic.
}
\usage{
findLargest(gN, testStat, data = "hgu133plus2")
}

\arguments{
  \item{gN}{A vector of probe identifiers for the chip. }
  \item{testStat}{A vector of test statistics, of the same length as
    \code{gN} with the per probe test statistics. }
  \item{data}{The character string identifying the chip. }
}
\details{
  All the probe identifiers, \code{gN}, are mapped to Entrez Gene IDs
  and the duplicates determined.  For any set of probes that map to the
  same Gene ID, the one with the largest test statistic is found. The
  return vector is the named vector of selected probe identifiers. The
  names are the Entrez Gene IDs.

  This could be extended in different ways, such as allowing the user to
  use a different selection criterion.  Also, matching on different
  identifiers seems like another alternative.
}
\value{
 A named vector of probe IDs. The names are Entrez Gene IDs.
}

\author{R. Gentleman }

\seealso{\code{\link{sapply}} }
\examples{
  library("hgu95av2")
  set.seed(124)
  gN = sample(ls(hgu95av2ENTREZID), 200)
  stats = rnorm(200)
  findLargest(gN, stats, "hgu95av2")

}
\keyword{manip}

