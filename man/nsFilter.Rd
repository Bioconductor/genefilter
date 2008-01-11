\name{nsFilter}

\alias{nsFilter}
\alias{varFilter}
\alias{featureFilter}
\alias{nsFilter,ExpressionSet-method}

\title{Non-Specific-ly Filter an ExpressionSet}
\description{
  This function identifies and removes probesets that are unlikely
  to be of use when modeling the data.  No phenotype variables are used
  in the filtering process, so the result can be used with any downstream
  analysis.
}
\usage{
nsFilter(eset, require.entrez = TRUE, require.GOBP = FALSE, 
    require.GOCC = FALSE, require.GOMF = FALSE,
    remove.dupEntrez = TRUE, var.func = IQR, var.cutoff = 0.5, 
    var.filter = TRUE, feature.exclude="^AFFX", ...)
varFilter(eset, var.func = IQR, var.cutoff = 0.5)
featureFilter(eset, require.entrez=TRUE,
    require.GOBP=FALSE, require.GOCC=FALSE,
    require.GOMF=FALSE, remove.dupEntrez=TRUE,
    feature.exclude="^AFFX")
}

\arguments{
  \item{eset}{an \code{ExpressionSet} object}
  \item{require.entrez}{If \code{TRUE}, require that all probe sets
      have an Entrez Gene ID annotation.  Probe sets without such an
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
  \item{var.func}{A \code{function} that will be used to assess the
      variance of a probe set across all samples.  This function
      should return a numeric vector of length one when given a
      numeric vector as input.  Probe sets with a \code{var.func}
      value less than \code{var.cutoff} will be removed. The default
      is \code{IQR}.}
  \item{var.cutoff}{A numeric value to use in filtering out probe sets
      with small variance across samples.  See the \code{var.func}
      argument and the details section below.}
  \item{var.filter}{A logical indicating whether or not to perform
      variance based filtering.  The default is \code{TRUE}.  A value of
      \code{FALSE} implies \code{remove.dupEntrez} is also forced to
      \code{FALSE}.}
  \item{feature.exclude}{A character vector of regular expressions.  Any
    probe sets identifiers (return value of \code{featureNames(eset)})
    that match one of the specified patterns will be filtered out.  The
    default value is intended to filter out Affymetrix quality control
    probe sets.}
  \item{...}{Unused, but available for specializing methods.}
}
\details{
  A first step in many microarray analysis procedures is to carry out
  non-specific filtering.  The goal is to remove uninteresting probe
  sets without regard to the phenotype data and reduce the number of
  probe sets that will be included in further analysis.

  \emph{Annotation Based Filtering} Arguments \code{require.entrez},
  \code{require.GOBP}, \code{require.GOCC}, and \code{require.GOMF}
  turn on a filter based on available annotation data.  The annotation
  package is determined by calling \code{annotation(eset)}.

  \emph{Variance Based Filtering} The \code{var.func} and
  \code{var.cutoff} arguments control the variance based filtering.  The
  intention is to remove probe sets with little variation across
  samples.  The default \code{var.func} is \code{IQR} and was selected
  because it is robust to outliers.  The deafult \code{var.cutoff} is
  \code{0.5} and is motivated by the common case where the platform is a
  genome-wide expression array and the rule of thumb that in any given
  tissue only 40\% of genes are expressed.

  Note that \code{IQR} is a reasonable choice only for the case of
  similar sample sizes. If one sample group is much bigger than the
  other(s), the values of the smaller group(s) will be discarded as
  outliers, and the \code{IQR} will be dominated by the large
  group. Hence the probe set is filtered out for low variance. One could
  consider using non-robust measures of variance (like \code{sd}) for
  these cases, or to look at the range of more extreme quantiles, where
  the quantile should be much smaller/larger than the ratios between
  sample sizes.

  \code{varFilter} does only variance filtering, and returns an
  \code{ExpressionSet}. \code{featureFilter} does only feature based
  filtering and returns an expression set.

}
\value{
 For \code{nsFilter} a list consisting of:
  \item{eset}{the filtered \code{ExpressionSet}}
  \item{filter.log}{a list giving details of how many probe sets where
    removed for each filtering step performed.}

  For both \code{varFilter} and \code{featureFilter} the filtered
  \code{ExpressionSet}.
}

\author{Seth Falcon}

\examples{
library("hgu95av2")
data(sample.ExpressionSet)
ans <- nsFilter(sample.ExpressionSet)
ans$eset
ans$filter.log

## skip variance-based filtering
ans <- nsFilter(sample.ExpressionSet, var.filter=FALSE)

a1 = varFilter(sample.ExpressionSet)
a2 = featureFilter(sample.ExpressionSet)
}


\keyword{ manip }

