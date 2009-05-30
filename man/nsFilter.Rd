\name{nsFilter}

\alias{nsFilter}
\alias{varFilter}
\alias{featureFilter}
\alias{nsFilter,ExpressionSet-method}

\title{Non-Specific Filtering of Features in an ExpressionSet}
\description{
  This function identifies and removes features that appear
  to be less informative. Use cases for this function are:
  variable selection for subsequent sample clustering or classification
  tasks;
  independent filtering of features used in subsequent hypothesis
  testing, with the aim of increasing the detection rate
  (please see Details).
}
\usage{
nsFilter(eset, require.entrez = TRUE, require.GOBP = FALSE, 
    require.GOCC = FALSE, require.GOMF = FALSE,
    remove.dupEntrez = TRUE, var.func = IQR, var.cutoff = 0.5, 
    var.filter = TRUE, filterByQuantile=TRUE,
    feature.exclude="^AFFX", ...)
varFilter(eset, var.func = IQR, var.cutoff = 0.5, filterByQuantile=TRUE)
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
      variance based filtering.  The default is \code{TRUE}.}  
\item{filterByQuantile}{Logical: whether the variance-filter cutoff threshold
  should be interpreted as a quantile. Defaults to \code{TRUE}; if set
  to  \code{FALSE} the cutoff value is used directly ``as is''.}
\item{feature.exclude}{A character vector of regular expressions.  Any
    probe sets identifiers (return value of \code{featureNames(eset)})
    that match one of the specified patterns will be filtered out.  The
    default value is intended to filter out Affymetrix quality control
    probe sets.}
  \item{...}{Unused, but available for specializing methods.}
}
\details{
  \emph{Marginal type I errors}:
  Independent filtering of features used in subsequent hypothesis
  testing can increase the detection rate at the same
  \emph{marginal} type I error, as detailed in the following.
  Call \eqn{U^I}{U^1} the stage 1 filter statistic,
  \eqn{U^{II}}{U^2} the stage 2 test statistic for differential expression.
  Sufficient conditions for marginal type-I error control are:
  \itemize{
    \item \eqn{U^I}{U^1} the overall (across all samples) variance or
    mean, and
    \eqn{U^{II}}{U^2} the t-statistic (or any other scale and location
    invariant statistic);

    \item \eqn{U^I}{U^1} the overall mean, and
    \eqn{U^{II}}{U^2} the moderated t-statistic (as in limma's \code{\link[limma:ebayes]{eBayes}}
    function);
    
    \item \eqn{U^I}{U^1} a sample-class label independent function
    (e.g. overall mean, median, variance, IQR), and
    \eqn{U^{II}}{U^2} the Wilcoxon rank sum statistic.
  }
  In each of these cases, the value of \eqn{U^I}{U^1} for the k-th
  feature must depend on the data for the k-th feature only, not on any
  other features.
  
  \emph{Experiment-wide type I error}:
  Marginal type-I error control provided by the conditions above
  is sufficient for control of the family wise error rate (FWER).
  Note, however, that common false discovery rate (FDR) methods depend 
  not only on the marginal behaviour of the test statistics under the
  null hypothesis, but also on their joint distribution.
  The joint distribution can be affected by filtering. The
  effect of this is negligible in many cases in practice, but this
  depends on the dataset and the filter used, and the assessment 
  is in the responsibility of the data analyst.
  For a more comprehensive discussion, please see the reference (Bourgon et al.).
  
  \emph{Annotation Based Filtering} Arguments \code{require.entrez},
  \code{require.GOBP}, \code{require.GOCC}, and \code{require.GOMF}
  turn on a filter based on available annotation data.  The annotation
  package is determined by calling \code{annotation(eset)}.

  \emph{Duplicate Probe Removal} If \code{remove.dupEntrez=TRUE},
  probes determined by your annotation to be pointing to the same gene
  will be compared, and only the probe with the highest \code{var.func} value
  will be retained.
  
  \emph{Variance Based Filtering} The \code{var.filter},
  \code{var.func}, \code{var.cutoff} and \code{varByQuantile} arguments
  control numerical cutoff-based filtering.  The intention is to remove
  uninformative probe sets, representing genes that were not expressed
  at all. Probes for which \code{var.func} returns \code{NA} are
  removed. The default \code{var.func} is \code{IQR}, which is defined as
  \code{rowQ(eset, ceiling(0.75 * ncol(eset))) - rowQ(eset, floor(0.25 * ncol(eset)))};
  this choice is motivated by the observation that unexpressed genes are
  detected most reliably through their low variability across samples.
  Additionally, \code{IQR} is robust to outliers (see note below). The
  default \code{var.cutoff} is \code{0.5} and is motivated by the rule of
  thumb that in many tissues only 40\% of genes are expressed. Of course,
  if you believe in a different approach to numerical filtering you can
  choose another function as \code{var.func}, or turn off numerical
  filtering by setting \code{var.filter=FALSE}.

  Note that by default the numerical-filter cutoff is interpreted
  as a quantile, so leaving the default values intact would filter out
  50\% of the genes remaining at this stage. If you prefer to set the
  cutoff at some absolute threshold, change the value of
  \code{varByQuantile} to \code{FALSE}, and modify \code{var.cutoff}
  accordingly.

  Note also that variance filtering is performed last, so that
  (if \code{varByQuantile=TRUE} and \code{remove.dupEntrez=TRUE}) the
  final number of genes does indeed exclude precisely the \code{var.cutoff} 
  fraction of unique genes remaining after all other filters were
  passed.
  
  The stand-alone function \code{varFilter} does only numerical filtering,
  and returns an \code{ExpressionSet}. \code{featureFilter} does only
  feature based filtering and duplicate removal, and returns an expression
  set as well. Duplicate removal is hard-coded to retain the highest-IQR
  probe for each gene.

}
\value{
  For \code{nsFilter} a list consisting of:
  \item{eset}{the filtered \code{ExpressionSet}}
  \item{filter.log}{a list giving details of how many probe sets where
    removed for each filtering step performed.}

  For both \code{varFilter} and \code{featureFilter} the filtered
  \code{ExpressionSet}.
}

\author{Seth Falcon (somewhat revised by Assaf Oron)}

\note{\code{IQR} is a reasonable variance-filter choice when the dataset
  is split into two roughly equal and relatively homogeneous phenotype
  groups. If your dataset has important groups smaller than 25\% of the
  overall sample size, or if you are interested in unusual
  individual-level patterns, then \code{IQR} may not be sensitive enough
  for your needs. In such cases, you should consider using less robust
  and more sensitive measures of variance (the simplest of which would
  be \code{sd}).}

\references{
  R. Bourgon, R. Gentleman, W. Huber,
  Independent filtering increases power for detecting differentially
  expressed genes, Technical Report.
}

\examples{
  library("hgu95av2.db")
  library("Biobase")
  data(sample.ExpressionSet)
  ans <- nsFilter(sample.ExpressionSet)
  ans$eset
  ans$filter.log

  ## skip variance-based filtering
  ans <- nsFilter(sample.ExpressionSet, var.filter=FALSE)

  a1 <- varFilter(sample.ExpressionSet)
  a2 <- featureFilter(sample.ExpressionSet)
}

\keyword{manip}
