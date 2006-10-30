\name{rowROC-class}
\docType{class}
\alias{rowROC-class}
\alias{pAUC}
\alias{AUC}
\alias{sens}
\alias{spec}
\alias{area}
\alias{pAUC,rowROC,numeric-method}
\alias{plot,rowROC,missing-method}
\alias{AUC,rowROC-method}
\alias{spec,rowROC-method}
\alias{sens,rowROC-method}
\alias{area,rowROC-method}
\alias{show,rowROC-method}
\alias{[,rowROC-method}


\title{Class "rowROC"}
\description{A class to model ROC curves and corresponding area under
  the curve as produced by rowpAUCs.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rowROC", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{data}:}{Object of class \code{"matrix"} The input data.}
    \item{\code{ranks}:}{Object of class \code{"matrix"} The ranked
      input data. }
    \item{\code{sens}:}{Object of class \code{"matrix"} Matrix of
      senitivity values for each gene at each cutpoint. }
    \item{\code{spec}:}{Object of class \code{"matrix"} Matrix of
      specificity values for each gene at each cutpoint.}
    \item{\code{pAUC}:}{Object of class \code{"numeric"} The partial
      area under the curve (integrated from 0 to \code{p}. }
    \item{\code{AUC}:}{Object of class \code{"numeric"} The total area
      under the curve. }
    \item{\code{factor}:}{Object of class \code{"factor"} The factor
      used for classification.}
    \item{\code{cutpoints}:}{Object of class \code{"matrix"} The values
      of the cutpoints at which specificity ans sensitivity was
calculated. (Note: the data is ranked prior to computation
      of ROC curves, the cutpoints map to the ranked data.}
    \item{\code{caseNames}:}{Object of class \code{"character"} The
      names of the two classification cases.}
    \item{\code{p}:}{Object of class \code{"numeric"} The limit to which
\code{pAUC} is integrated. }
  }
}
\section{Methods}{
  \describe{
    \item{show \code{signature(object="rowROC")}}{Print nice info
      about the object.}
     \item{[ \code{signature(x="rowROC", j="missing")}}{Subset the
       object according to rows/genes.}
     \item{plot \code{signature(x="rowROC", y="missing")}}{Plot the ROC
       curve of the first row of the object along with the \code{pAUC}.
     To plot the curve for a specific row/gene subsetting should be done
     first (i.e. \code{plot(rowROC[1])}.}
    \item{pAUC \code{signature(object="rowROC", p="numeric")}}{Integrate
      area under the curve from \code{0} to \code{p}. This method
      returns a new \code{rowROC} object.}
    \item{AUC \code{signature(object="rowROC")}}{Integrate
      total area under the curve. This method returns a new
      \code{rowROC} object.}
    \item{sens \code{signature(object="rowROC")}}{Accessor method for
      sensitivity slot.}
    \item{spec \code{signature(object="rowROC")}}{Accessor method for
      specificity slot.}
    \item{area \code{signature(object="rowROC", total="logical")}}{Accessor method for
       pAUC slot.}
  }
}
\references{Pepe MS, Longton G, Anderson GL,
Schummer M.: Selecting
differentially expressed genes from microarray
experiments. \emph{Biometrics. 2003 Mar;59(1):133-42.}}
\author{Florian Hahne <f.hahne@dkfz.de>}


\seealso{
  \code{\link[genefilter:rowpAUCs]{rowpAUCs}}
}
\examples{
 data(sample.ExpressionSet)
roc <- rowpAUCs(sample.ExpressionSet, "sex", p=0.5)
roc
area(roc[1:3])

if(interactive()) {
par(ask=TRUE)
plot(roc)
plot(1-spec(roc[1]), sens(roc[2]))
par(ask=FALSE)
}

pAUC(roc, 0.1)
roc
}
\keyword{classes}
