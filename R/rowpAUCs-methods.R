## ==========================================================================
## core rowpAUCs method for objects of class matrix
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("rowpAUCs", signature(x="matrix", fac="factor"),
          function(x, fac, p=0.1, flip=TRUE, caseNames=c("1", "2")){  
            ##check argument 'p'
            if(!is.numeric(p) || length(p)>1)
              stop("'p' must be numeric of length 1")
            ## check argument 'fac'
            f <- genefilter:::checkfac(fac)
            if(f$nrgrp != 2 || length(f$fac) != ncol(x) ||
               length(unique(f$fac)) !=2 )
                stop("'fac' must be factor with 2 levels and length 'ncol(x)'")
            ## check argument 'flip'
            if(length(flip)!=1 || !(is.logical(flip)))
                stop("'flip' must be logical scalar")
            flip <- as.integer(flip)
            ## compute cutpoints
            cutpts <- matrix((0:ncol(x))+0.5, ncol=ncol(x)+1, 
                             nrow=nrow(x), byrow=TRUE,
                             dimnames=list(rownames(x), NULL))
            
            ## rank data
            xr <- t(apply(x, 1, rank))
            mode(xr) <- "numeric"
            ## call C function and return object of class 'rowROC'
            res <- .Call("ROCpAUC", xr, cutpts, as.integer(f$fac), p,
                         PACKAGE="genefilter", flip)
            sens <- res$sens
            spec <- res$spec
            rownames(sens) <- rownames(spec) <- rownames(x)
            pAUC <- res$pAUC
            AUC <- res$AUC
            names(AUC) <- names(pAUC) <- rownames(x)
            object <- new("rowROC", data=x, sens=sens,
                          spec=spec, pAUC=pAUC, AUC=AUC,
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
          function(x, fac, p=0.1, flip=TRUE, caseNames=c("1", "2")){
            cutpts <- matrix((0:ncol(x))+0.5, ncol=ncol(x)+1, 
                             nrow=nrow(x), byrow=TRUE)
            rowpAUCs(x=x, fac=factor(fac), p=p, flip=flip,
                   caseNames=caseNames) 
          })
## ==========================================================================


## ==========================================================================
## rowpAUCs method with signature x=ExpressionSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
setMethod("rowpAUCs", signature(x="ExpressionSet"),        
          function(x, fac, p=0.1, flip=TRUE, caseNames=c("1", "2")){
            rowpAUCs(x=exprs(x), fac=fac, p=p, flip=flip,
                     caseNames=caseNames) 
          })
## ==========================================================================


## ==========================================================================
## rowpAUCs method with signature x=ExpressionSet fac=character
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
setMethod("rowpAUCs", signature(x="ExpressionSet", fac="character"),        
          function(x, fac, p=0.1, flip=TRUE, caseNames=c("1", "2")){
            if (length(fac) == 1){
                if(!fac %in% colnames(pData(x)))
                    stop("fac must be length 1 character indicating a ",
                         "covariate in the phenoData slot of the expressionSet")
                cn <- as.character(levels(pData(x)[[fac]]))
                fac = factor(as.integer(factor(pData(x)[[fac]]))-1)
                rowpAUCs(x=exprs(x), fac=fac, p=p, flip=flip, caseNames=cn)
            }else{
                rowpAUCs(x=x, fac=as.factor(fac), p=p, flip=flip,
                         caseNames=caseNames)
            }
          })
## ==========================================================================
