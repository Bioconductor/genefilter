\name{rowpAUCs-methods}
\docType{methods}
\alias{rowpAUCs-methods}
\alias{rowpAUCs}
\alias{rowpAUCs,matrix,factor-method}
\alias{rowpAUCs,matrix,numeric-method}
\alias{rowpAUCs,ExpressionSet,ANY-method}
\alias{rowpAUCs,ExpressionSet,character-method}
\title{Rowwise ROC and pAUC computation}
\description{Methods for fast rowwise computation of ROC curves and
  (partial) area under the curve (pAUC) using the simple classification
  rule \code{x > theta}, where \code{theta} is a value in the range of
  \code{x}
}
\usage{
rowpAUCs(x, fac, p=0.1, flip=TRUE, caseNames=c("1", "2"))
}
\arguments{
  \item{x}{\code{ExpressionSet} or numeric \code{matrix}. The
    \code{matrix} must not contain \code{NA} values.}

  \item{fac}{A \code{factor} or \code{numeric} or \code{character} that can
    be coerced to a \code{factor}. If \code{x} is an \code{ExpressionSet},
    this may also be a character \code{vector} of length 1 with the name of
    a covariate variable in \code{x}. \code{fac} must have exactly 2 levels.
    For better control over the classification, use integer values in 0 and 1,
    where 1 indicates the "Disease" class in the sense of the Pepe et al paper
    (see below).}
 
  \item{p}{Numeric \code{vector} of length 1. Limit in (0,1) to integrate pAUC
    to.}
  
  \item{flip}{Logical. If \code{TRUE}, both classification rules \code{x
  > theta} and \code{x < theta} are tested and the (partial) area under
  the curve of the better one of the two is returned. This is
  appropriate for the cases in which the classification is not
  necessarily linked to higher expression values, but instead it is
  symmetric and one would assume both over- and under-expressed genes for
  both classes. You can set \code{flip} to \code{FALSE} if you only want
  to screen for genes which discriminate Disease from Control with the
  \code{x > theta} rule.}

  \item{caseNames}{The class names that are used when plotting the
    data. If \code{fac} is the name of the covariate variable in the
    \code{ExpressionSet} the function will use its levels as
    \code{caseNames}.}
}
\details{
  Rowwise calculation of Receiver Operating Characteristic (ROC) curves
  and the corresponding partial area under the curve (pAUC) for a given
  data matrix or \code{ExpressionSet}. The function is implemented in C
  and thus reasonably fast and memory efficient. Cutpoints (\code{theta}
  are calculated before the first, in between and after the last data
  value. By default, both classification rules \code{x > theta} and
  \code{x < theta} are tested and the (partial) area under the curve of
  the better one of the two is returned. This is only valid for
  symmetric cases, where the classification is independent of the
  magnitude of \code{x} (e.g., both over- and under-expression of
  different genes in the same class).  For unsymmetric cases in which
  you expect x to be consistently higher/lower in of of the two classes
  (e.g. presence or absence of a single biomarker) set \code{flip=FALSE}
  or use the functionality provided in the \code{ROC} package. For
  better control over the classification (i.e., the choice of "Disease"
  and "Control" class in the sense of the Pepe et al paper), argument
  \code{fac} can be an integer in \code{[0,1]} where 1 indicates
  "Disease" and 0 indicates "Control". 
}
\section{Methods}{
  \describe{
    Methods exist for \code{rowPAUCs}:
    \item{rowPAUCs}{\code{signature(x="matrix", fac="factor")}}
    \item{rowPAUCs}{\code{signature(x="matrix", fac="numeric")}}
    \item{rowPAUCs}{\code{signature(x="ExpressionSet")}}
    \item{rowPAUCs}{\code{signature(x="ExpressionSet", fac="character")}}
  }
}
\value{
  An object of class \code{\link[genefilter:rowROC]{rowROC}} with the
  calculated specificities and sensitivities for each row and the
  corresponding pAUCs and AUCs values. See
  \code{\link[genefilter:rowROC]{rowROC}} for details.
}
\references{Pepe MS, Longton G, Anderson GL,
Schummer M.: Selecting
    differentially expressed genes from microarray experiments. \emph{Biometrics. 2003 Mar;59(1):133-42.}}
\author{Florian Hahne <fhahne@fhcrc.org>}
\seealso{\code{\link[ROC:rocdemo.sca]{rocdemo.sca},
    \link[ROC:AUC]{pAUC}, \link[genefilter:rowROC]{rowROC}}}
\examples{
data(sample.ExpressionSet)

r1 = rowttests(sample.ExpressionSet, "sex")
r2 = rowpAUCs(sample.ExpressionSet, "sex", p=0.1)

plot(area(r2, total=TRUE), r1$statistic, pch=16)
sel <- which(area(r2, total=TRUE) > 0.7)
plot(r2[sel])

## this compares performance and output of rowpAUCs to function pAUC in
## package ROC 
if(require(ROC)){
  ## performance
  myRule = function(x)
    pAUC(rocdemo.sca(truth = as.integer(sample.ExpressionSet$sex)-1 ,
         data = x, rule = dxrule.sca), t0 = 0.1)
  nGenes = 200
  cat("computation time for ", nGenes, "genes:\n")
  cat("function pAUC: ")
  print(system.time(r3 <- esApply(sample.ExpressionSet[1:nGenes, ], 1, myRule)))
  cat("function rowpAUCs: ")
  print(system.time(r2 <- rowpAUCs(sample.ExpressionSet[1:nGenes, ],
  "sex", p=1)))

  ## compare output
  myRule2 = function(x)
   pAUC(rocdemo.sca(truth = as.integer(sample.ExpressionSet$sex)-1 ,
                    data = x, rule = dxrule.sca), t0 = 1)
  r4 <-  esApply(sample.ExpressionSet[1:nGenes, ], 1, myRule2)
  plot(r4,area(r2), xlab="function pAUC", ylab="function rowpAUCs",
  main="pAUCs")

  plot(r4, area(rowpAUCs(sample.ExpressionSet[1:nGenes, ],
  "sex", p=1, flip=FALSE)), xlab="function pAUC", ylab="function rowpAUCs",
  main="pAUCs")

  r4[r4<0.5] <- 1-r4[r4<0.5]
  plot(r4, area(r2), xlab="function pAUC", ylab="function rowpAUCs",
  main="pAUCs")
 }
}
\keyword{math}
