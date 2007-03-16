\name{nsFilter}

\alias{nsFilter}
\alias{nsFilter,ExpressionSet-method}

\title{Non-Specific-ly Filter an ExpressionSet}
\description{
  This function removes unwanted probe sets from an
  \code{ExpressionSet} without using phenotype variables in the
  filtering process.  Hence the filter is non-specific with respect to
  the phenotypes in the data.
}
\usage{
nsFilter(eset, require.entrez = TRUE, require.symbol = TRUE,
require.GOBP = FALSE, require.GOCC = FALSE, require.GOMF = FALSE,
remove.dupEntrez = TRUE, var.func = IQR, var.cutoff = 0.5, var.filter = TRUE)
}

\arguments{
  \item{eset}{an \code{ExpressionSet} object}
  \item{require.entrez}{If \code{TRUE}, require that all probe sets
      have an Entrez Gene ID annotation.  Probe sets without such an
      annotation will be filtered out.}
  \item{require.symbol}{If \code{TRUE}, require that all probe sets
      have a gene symbol annotation.  Probe sets without such an
      annotation will be filtered out.}
  \item{require.GOBP}{If \code{TRUE}, require that all probe sets have
    an annotation to at least one GO ID in the BP ontology.  Probe
    sets without such an annotation will be filtered out.}
  \item{require.GOCC}{If \code{TRUE}, require that all probe sets have
    an annotation to at least one GO ID in the CC ontology.  Probe
    sets without such an annotation will be filtered out.}
  \item{require.GOMF}{If \code{TRUE}, require that all probe sets have
    an annotation to at least one GO ID in the MF ontology.  Probe
    sets without such an annotation will be filtered out.}
  \item{remove.dupEntrez}{If \code{TRUE} and there are multiple probe
      sets mapping to the same Entrez Gene ID, then the probe set with
      the largest value of \code{var.func} will be retained and the
      others removed.}
  \item{var.func}{a \code{function} that will be used to assess the
      variance of a probe set across all samples.  This function
      should return a numeric vector of length one when given a
      numeric vector as input.  Probe sets with a \code{var.func}
      value less than \code{var.cutoff} will be removed. The default
      is \code{IQR}.}
  \item{var.cutoff}{a numeric value to use in filtering out probe sets
      with small variance across samples.  See the \code{var.func}
      argument and the details section below.}
  \item{var.filter}{a logical indicating whether or not to perform
      variance based filtering.  The default is \code{TRUE}.}
}
\details{
  A first step in many microarray analysis procedures is to carry out
  non-specific filtering.  The goal is to remove uninteresting probe
  sets without regard to the phenotype data and reduce the number of
  probe sets that will be included in further analysis.

  \emph{Annotation Based Filtering}
  Arguments \code{require.entrez}, \code{require.symbol},
  \code{require.GOBP}, \code{require.GOCC}, and \code{require.GOMF}
  turn on a filter based on available annotation data.  The annotation
  package is determined by calling \code{annotation(eset)}.

  \emph{Variance Based Filtering}
  The \code{var.func} and \code{var.cutoff} arguments control the
  variance based filtering.  The intention is to remove probe sets
  with little variation across samples.  The default \code{var.func}
  is \code{IQR} and was selected because it is robust to outliers.
  The deafult \code{var.cutoff} is \code{0.5} and is motivated by the
  common case where the platform is a genome-wide expression array and
  the rule of thumb that in any given tissue only 40\% of genes are
  expressed.

}
\value{
  A list consisting of:
  \item{eset}{the filtered \code{ExpressionSet}}
  \item{filter.log}{a list giving details of how many probe sets where
    removed for each filtering step performed.}
}

\author{Seth Falcon}

\examples{
library("hgu95av2")
data(sample.ExpressionSet)
ans <- nsFilter(sample.ExpressionSet)
ans$eset
ans$filter.log
}


\keyword{ manip }

