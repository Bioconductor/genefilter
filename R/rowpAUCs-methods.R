## ==========================================================================
## core rowpAUCs method for objects of class matrix
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("rowpAUCs", signature(x="matrix", fac="factor"),
          function(x, fac, p=0.1, jitter=TRUE, caseNames=c("1", "2")){  
            ##check argument 'p'
            if(!is.numeric(p) || length(p)>1)
              stop("'p' must be numeric of length 1")
            ## check argument 'fac'
            f = genefilter:::checkfac(fac)
            if(f$nrgrp != 2 || length(f$fac) != ncol(x) ||
               length(unique(f$fac)) !=2 )
              stop("'fac' must be factor with 2 levels and length 'ncol(x)'")
            ## compute cutpoints
            cutpts <- matrix((0:ncol(x))+0.5, ncol=ncol(x)+1, 
                             nrow=nrow(x), byrow=TRUE)
            ## rank data
            xr <- t(apply(x, 1, rank))
            mode(xr) <- "numeric"
            ## add jitter 
            if(jitter)
              xr <- jitter(xr, amount=0.1)
            ## call C function and return object of class 'rowROC'
            res <- .Call("ROCpAUC", xr, cutpts, as.integer(f$fac), p,
                         PACKAGE="genefilter")
            object <- new("rowROC", data=x, sens=res$sens,
                          spec=res$spec, pAUC=res$pAUC, AUC=res$AUC,
                          factor=factor(f$fac), p=p, ranks=xr,
                          caseNames=as.character(caseNames),
                          cutpoints=cutpts)
            return(object)
          })
## ==========================================================================


## ==========================================================================
## rowpAUCs method with signature x=matrix, fac=numeric
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("rowpAUCs", signature(x="matrix", fac="numeric"),
          function(x, fac, p=0.1, jitter=FALSE, caseNames=c("1", "2")){
            cutpts <- matrix((0:ncol(x))+0.5, ncol=ncol(x)+1, 
                             nrow=nrow(x), byrow=TRUE)
            rowpAUCs(x=x, fac=factor(fac), p=p, jitter=jitter,
                   caseNames=caseNames) 
          })
## ==========================================================================


## ==========================================================================
## rowpAUCs method with signature x=exprSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("rowpAUCs", signature(x="exprSet"),        
          function(x, fac, p=0.1, jitter=FALSE, caseNames=c("1", "2")){
            .Deprecated(msg=Biobase:::EXPRSET_DEPR_MSG)
            rowpAUCs(x=exprs(x), fac=fac, p=p, jitter=jitter,
                     caseNames=caseNames) 
          })

## ==========================================================================


## ==========================================================================
## rowpAUCs method with signature x=exprSet fac=character
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
setMethod("rowpAUCs", signature(x="exprSet", fac="character"),        
          function(x, fac, p=0.1, jitter=FALSE, caseNames=c("1", "2")){
            .Deprecated(msg=Biobase:::EXPRSET_DEPR_MSG)
            if (length(fac) != 1)
               stop("fac must be length 1 character or a factor")
            cn <- as.character(levels(pData(x)[[fac]]))
            fac = factor(as.integer(factor(pData(x)[[fac]]))-1)
            rowpAUCs(x=exprs(x), fac=fac, p=p, jitter=jitter, caseNames=cn)   
          })
## ==========================================================================


## ==========================================================================
## rowpAUCs method with signature x=ExpressionSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
setMethod("rowpAUCs", signature(x="ExpressionSet"),        
          function(x, fac, p=0.1, jitter=FALSE, caseNames=c("1", "2")){
            rowpAUCs(x=exprs(x), fac=fac, p=p, jitter=jitter,
                     caseNames=caseNames) 
          })
## ==========================================================================


## ==========================================================================
## rowpAUCs method with signature x=ExpressionSet fac=character
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
setMethod("rowpAUCs", signature(x="ExpressionSet", fac="character"),        
          function(x, fac, p=0.1, jitter=FALSE, caseNames=c("1", "2")){
            if (length(fac) != 1)
               stop("fac must be length 1 character or a factor")
            cn <- as.character(levels(pData(x)[[fac]]))
            fac = factor(as.integer(factor(pData(x)[[fac]]))-1)
            rowpAUCs(x=exprs(x), fac=fac, p=p, jitter=jitter, caseNames=cn) 
          })
## ==========================================================================
