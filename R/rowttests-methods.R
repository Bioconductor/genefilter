setGeneric("rowttests", function(x, fac, tstatOnly=FALSE)
           standardGeneric("rowttests"))

## Core method for matrix

setMethod("rowttests", signature(x="matrix", fac="factor"),
          function(x, fac, tstatOnly=FALSE) {
              if (!missing(tstatOnly)
                  && (!is.logical(tstatOnly) || is.na(tstatOnly)))
                stop(sQuote("tstatOnly"), " must be TRUE or FALSE.")
              f <- checkfac(fac)
              if (f$nrgrp > 2)
                stop("Number of groups must be <= 2 for 'rowttests'.")
              res <- .Call("rowcolttests", x, f$fac, f$nrgrp,
                           as.integer(0), PACKAGE="genefilter")
              if (!tstatOnly)
                res$p.value <- 2 * pt(abs(res$statistic), res$df,
                                      lower.tail=FALSE)
              return(res)
          })


setMethod("rowttests", signature(x="matrix", fac="missing"),
          function(x, fac, tstatOnly=FALSE) {
              rowttests(x, fac=factor(integer(ncol(x))), tstatOnly)
          })


## exprSet methods

setMethod("rowttests", signature(x="exprSet", fac="missing"),
          function(x, fac, tstatOnly=FALSE) {
              .Deprecated(msg=Biobase:::EXPRSET_DEPR_MSG)           
              x <- exprs(x)
              fac <- integer(ncol(x))
              rowttests(x, fac=fac, tstatOnly=tstatOnly)
          })


setMethod("rowttests", signature(x="exprSet", fac="factor"),
          function(x, fac, tstatOnly=FALSE) {
              .Deprecated(msg=Biobase:::EXPRSET_DEPR_MSG)
              x <- exprs(x)
              rowttests(x, fac=fac, tstatOnly=tstatOnly)
          })


setMethod("rowttests", signature(x="exprSet", fac="character"),
          function(x, fac, tstatOnly=FALSE) {
              .Deprecated(msg=Biobase:::EXPRSET_DEPR_MSG)
              if (length(fac) != 1)
                stop("fac must be length 1 character or a factor")
              fac <- factor(pData(x)[[fac]])
              rowttests(exprs(x), fac=fac, tstatOnly=tstatOnly)
          })


## ExpressionSet methods

setMethod("rowttests", signature(x="ExpressionSet", fac="missing"),
          function(x, fac, tstatOnly=FALSE) {
              x <- exprs(x)
              fac <- integer(ncol(x))
              rowttests(x, fac=fac, tstatOnly=tstatOnly)
          })


setMethod("rowttests", signature(x="ExpressionSet", fac="factor"),
          function(x, fac, tstatOnly=FALSE) {
              x <- exprs(x)
              rowttests(x, fac=fac, tstatOnly=tstatOnly)
          })


setMethod("rowttests", signature(x="ExpressionSet", fac="character"),
          function(x, fac, tstatOnly=FALSE) {
              if (length(fac) != 1)
                stop("fac must be length 1 character or a factor")
              fac <- factor(pData(x)[[fac]])
              rowttests(exprs(x), fac=fac, tstatOnly=tstatOnly)
          })

