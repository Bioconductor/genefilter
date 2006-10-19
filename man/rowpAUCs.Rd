\name{rowpAUCs}
\alias{rowpAUCs}

\title{Rowwise ROC and pAUC}
\description{Fast rowwise calculation of ROC curves and pAUC.}

\usage{
rowpAUCs(x, fac, cutpts, p=0.1, common=FALSE)
}
\arguments{
  \item{x}{\code{exprSet} or numeric matrix. The matrix must not
    contain \code{NA} values.}
  \item{fac}{Factor; if \code{x} is an \code{exprSet}, this may also be
    a character vector of length 1 with the name of a covariate
    variable in \code{x}. \code{fac} must have exactly 2 levels.}
  \item{cutpts}{Matrix with same number of rows as \code{x} or vector
    of length greater than 1. It specifies the thresholds for the
    calculation of the ROC curves. Smaller numbers may further speed up
    computations. If missing, the ROC curves are calculated between data
    points.}
  \item{p}{Numeric vector of length 1. Limit in (0,1) to integrate pAUC
    to.}
  \item{common}{Calculate common cutpoints for the whole dataset. If
    \code{FALSE} cutpoints are calculated for each row individually.}
}
\details{
  Rowwise calculation of Receiver Operating
  Characteristic (ROC) curves and the corresponding partial area under
  the curve (pAUC) for a given data matrix or \code{exprSet}. The
  function is implemented in C and thus reasonably fast and memory
  efficient.
}
\value{
  A list with the calculated specificities and sensitivities for each
  row as matrices, and the corresponding pAUCs.  
}
\references{Pepe MS, Longton G, Anderson GL, Schummer M.:  Selecting
    differentially expressed genes from microarray experiments. \emph{Biometrics. 2003 Mar;59(1):133-42.}}
\author{Florian Hahne <f.hahne@dkfz.de>}
\seealso{\code{\link[ROC:rocdemo.sca]{rocdemo.sca}, \link[ROC:AUC]{pAUC}}}
\examples{

   data(sample.exprSet.1)

   r1 = rowttests(sample.exprSet.1, "cov2")
   r2 = rowpAUCs(sample.exprSet.1, "cov2")

 if(interactive()) {
   plot(r2$pAUC, r1$statistic, pch=16)
   plot(1-r2$spec[1,], r2$sens[1,], pch=16)
 }

 ## this compares performance and output of rowpAUCs to function pAUC in
 ## package ROC 
 if(require(ROC)){

  myRule = function(x) {
    pAUC(rocdemo.sca(truth = pData(sample.exprSet.1)$cov2-1, data = x, rule = dxrule.sca), t0 = 0.1)
  }
  cat("computation time:\n")
  cat("function pAUC: ")
  print(system.time(r3 <- esApply(sample.exprSet.1[1:100, ], 1, myRule)))
  cat("function rowpAUCs: ")
  print(system.time(r2 <- rowpAUCs(sample.exprSet.1[1:100, ], "cov2")))
  plot(r3, r2$pAUC, xlab="function pAUC", ylab="function rowpAUCs", main="pAUCs")
 }
}


}
\keyword{math}
