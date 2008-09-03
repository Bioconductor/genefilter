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
                   ylab="sensitivity",
                   main=NULL,
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
              nn <-  names(area(x)[1])
              if(is.null(main))
                  main <- paste("ROC-Curve", ifelse(length(nn),
                                                    paste("(", nn, ")", sep=""), ""))
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
              atext <- paste("AUC:  ", signif(x@AUC[1],3))
              tw <- strwidth(atext)
              w <- diff(par("usr")[1:2])
              cex <- min(1, (w/2+w/10)/tw)
              th <- strheight(atext, cex=cex)*1.1
              if(x@p<1){
                  ptext <- paste("pAUC: ", signif(x@pAUC[1],3), " (p=", x@p,  ")",
                                 sep="")
                  tw <- max(tw, strwidth(ptext))
                  cex <- min(1, (w/2+w/10)/tw)
                  abline(v=x@p, col="darkblue", lty=2) 
                  text(x=1-tw*cex*1.1, y=0.02+th*cex, atext, pos=4, cex=cex)
                  text(x=1-tw*cex*1.1, y=0.02, ptext, pos=4, cex=cex)
              }else{
                  text(x=1-tw*cex*1.1, y=0.02, atext, pos=4, cex=cex)
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
              names(res$pAUC) <-  names(res$AUC) <- names(object@AUC)
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
              object@pAUC <- object@AUC
              object@p <- 1
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
