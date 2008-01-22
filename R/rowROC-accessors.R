## ==========================================================================
## show method for objects of class rowROC
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show", signature(object="rowROC"), function(object){
  cat("matrix of ROC curves for", nrow(object@data), "genes/rows",
      "with", max(0,ncol(object@cutpoints)), "cutpoints\n")
  cat("  size of class ", object@caseNames[1] ,": ",
      sum(object@factor==levels(object@factor)[1]), "\n", sep="")
  cat("  size of class ", object@caseNames[2] ,": ",
      sum(object@factor==levels(object@factor)[2]), "\n", sep="")
  cat("partial areas under curve calculated for p=",
      object@p, "\n", sep="")
})
## ==========================================================================


## ==========================================================================
## subsetting method for objects of class rowROC
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("[",
  signature="rowROC",
  definition=function(x, i, j="missing", drop="missing") {
    x@sens <- x@sens[i,,drop=FALSE]
    x@spec <- x@spec[i,,drop=FALSE]
    x@pAUC <- x@pAUC[i]
    x@AUC <- x@AUC[i]
    x@data <- x@data[i,,drop=FALSE]
    x@cutpoints <- x@cutpoints[i,,drop=FALSE]
    x@ranks <- x@ranks[i,,drop=FALSE]
    return(x)
           },
   valueClass="rowROC")
## ==========================================================================


## ==========================================================================
## plot method for objects of class rowROC
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("plot", signature(x="rowROC", y="missing"),
          function(x, pch=20, cex=0.7, xlab="1 - specificity",
                   ylab="sensitivity", main="ROC-Curve",
                   sub=paste("class ", x@caseNames[1], " (",
                     sum(x@factor==levels(x@factor)[1]),
                     " cases) | class ", x@caseNames[2], " (",
                     sum(x@factor==levels(x@factor)[2]), " cases)", sep=""),
                   ...){
            sx <- sort(1-x@spec[1,])
            sy <- sort(x@sens[1,])
            spx <- c(sx[sx<=x@p & sy>0],x@p)
            spy <- sy[sx<=x@p & sy>0]
            if(!length(spy)){
              spy <- 0
              spx <- c(0,spx)
            }
            spy <- c(spy, max(spy))
            len <- length(sx)
            plot(sx, sy, pch=pch, cex=cex, xlab=xlab,
                 ylab=ylab, main=main, sub=sub, ...)
            if(mean(x@data)==1 || all(sx==sy))
              polygon(c(0,1,1), c(0,0,1), col="#ececec", lty=0)
            else{
              rect(spx[-1], 0, spx[-1] - diff(spx),spy[-1],
                   col="#ececec", lty=0)       
              lines(sx, sy, type="s")
            }
            points(sx, sy, pch=pch, cex=cex, ...)
            lines(0:1, 0:1, lty=3, col="darkgray")
            abline(v=x@p, col="darkblue", lty=2)
            if(x@p<1){
              text(x=0.70, y=0.06, paste("AUC:  ", signif(x@AUC[1],3)),pos=4)
              text(x=0.70, y=0.02, paste("pAUC: ", signif(x@pAUC[1],3),
                          " (p=", x@p,  ")", sep=""), pos=4)
            }else{
              text(x=0.75, y=0.02, paste("AUC:  ", signif(x@AUC[1],3)),pos=4)
            }
          })
## ==========================================================================


## ==========================================================================
## pAUC method for objects of class rowROC
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("pAUC", signature(object="rowROC", p="numeric"),
          function(object, p, flip=TRUE){
              if(length(flip)!=1 || !(is.logical(flip)))
                  stop("'flip' must be logical scalar")
              flip <- as.integer(flip)
              res <- .Call("pAUC", object@spec, object@sens, p, flip)
              object@pAUC <- res$pAUC
              object@AUC <- res$AUC
              object@p <- p
              return(object)
          })
## ==========================================================================


## ==========================================================================
## AUC method for objects of class rowROC
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("AUC", signature(object="rowROC"),          
          function(object){
          object <- pAUC(object, p=1)
          return(object)
        })
## ==========================================================================


## ==========================================================================
## accessor method to slot 'sens' for objects of class rowROC
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("sens", signature(object="rowROC"),          
          function(object)
          return(object@sens)
          )
## ==========================================================================


## ==========================================================================
## accessor method to slot 'spec' for objects of class rowROC
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("spec", signature(object="rowROC"),          
          function(object)
          return(object@spec)
          )
## ==========================================================================


## ==========================================================================
## accessor method to slots 'AUC' or 'pAUC' for objects of class rowROC
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("area", signature(object="rowROC"),          
          function(object, total=FALSE){
            if(total)
              return(object@AUC)
            else
              return(object@pAUC)
          })
## ==========================================================================
