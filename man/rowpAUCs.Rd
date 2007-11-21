\name{rowpAUCs-methods}
\docType{methods}
\alias{rowpAUCs-methods}
\alias{rowpAUCs}
\alias{rowpAUCs,matrix,factor-method}
\alias{rowpAUCs,matrix,numeric-method}
\alias{rowpAUCs,exprSet,ANY-method}
\alias{rowpAUCs,exprSet,character-method}
\alias{rowpAUCs,ExpressionSet,ANY-method}
\alias{rowpAUCs,ExpressionSet,character-method}

\title{Rowwise ROC and pAUC computation}
\description{Methods for fast rowwise computation of ROC curves and
  (partial) area under the curve (pAUC).}

\usage{
rowpAUCs(x, fac, p=0.1, jitter=FALSE, caseNames=c("1", "2"))
}
\arguments{
  \item{x}{\code{ExpressionSet}, \code{exprSet} or numeric matrix. The
matrix must not contain \code{NA} values. Note that the use of
\code{exprSet} is deprecated.}
  \item{fac}{Factor; if \code{x} is an \code{ExpressionSet} or
\code{exprSet}, this may also be a character vector of length 1 with the
name of a covariate variable in \code{x}. \code{fac} must have exactly 2
levels.}
 
  \item{p}{Numeric vector of length 1. Limit in (0,1) to integrate pAUC
    to.}
  \item{jitter}{Add small amount of jitter to the data. This can improve
    the look of the ROC curve when most of the values in \code{x} are
    actually the same since cutpoints are calculated between data
    values.}
  \item{caseNames}{caseNames}{ The
      names of the two classification cases.}
}
\details{
  Rowwise calculation of Receiver Operating
  Characteristic (ROC) curves and the corresponding partial area under
  the curve (pAUC) for a given data matrix or \code{ExpressionSet}. The
  function is implemented in C and thus reasonably fast and memory
  efficient. Cutpoints are calculated before the first, in between and
  after the last data value. Note: If the total area under the curve is
  smaller 0.5 it is flipped for plotting and computation of \code{pAUC}.
}
\section{Methods}{
  \describe{
    Methods exist for \code{rowPAUCs}:
    \item{rowPAUCs}{\code{signature(x="matrix", fac="factor")}}
    \item{rowPAUCs}{\code{signature(x="matrix", fac="numeric")}}
    \item{rowPAUCs}{\code{signature(x="exprSet")}}
    \item{rowPAUCs}{\code{signature(x="exprSet", fac="character")}}
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
\author{Florian Hahne <f.hahne@dkfz.de>}
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
  plot(r4,area(r2), xlab="function pAUC", ylab="function rowpAUCs", main="pAUCs")
  r4[r4<0.5] <- 1-r4[r4<0.5]
  plot(r4, area(r2), xlab="function pAUC", ylab="function rowpAUCs", main="pAUCs")
 }
}


\keyword{math}
