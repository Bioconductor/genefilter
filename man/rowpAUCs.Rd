\name{rowpAUCs}
\alias{rowpAUCs}

\title{Rowwise ROC and pAUC}
\description{Fast rowwise calculation of ROC curves and pAUC.}

\usage{
rowpAUCs(x, fac, cutpts, p=0.1)
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
  \item{p}{Numeric vector of length 1. Limit in (0,1) to integrate pAUC to.}
}
\details{
  Rowwise calculation of Receiver Operating
  Characteristic (ROC) curves and the corresponding partial areae under
  the curve (pAUC) for a given data matrix or \code{exprSet}. The
  function is implemented in C and thus reasonably fast and memory
  efficient.
  
  The definition of the pAUC uses a naive trapezoidal rule
  which, although less accurate, is considerably faster than real integration. 
}
\value{
  A list with the calculated specificities and sensitivities for each
  row as matrices, and the corresponding pAUCs.  
}
\references{Pepe MS, Longton G, Anderson GL, Schummer M.,  Selecting
    differentially expressed genes from microarray experiments. \emph{Biometrics. 2003 Mar;59(1):133-42.}}
\author{Florian Hahne <f.hahne@dkfz.de>}
\seealso{\code{\link[ROC:rocdemo.sca]{rocdemo.sca}, \link[ROC:pAUC]{pAUC}}}
\examples{
  
   x  = matrix(runif(9700), ncol=97)
   f2 = floor(runif(ncol(x))*2)

   r1 = rowttests(x, f2)
   r2 = rowpAUCs(x, f2)

 if(interactive())
   plot(r2$pAUC, r1$statistic, pch=16)
   plot(1-r2$spec[1,], r2$sens[1,], pch=16)   
}


}
\keyword{math}
