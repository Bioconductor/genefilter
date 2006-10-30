## Classes for package genefilter

## ==========================================================================
## class rowROC: objects model result of call to function rowpAUCs,
##               pAUC or AUC
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setClass("rowROC",
         representation(data = "matrix",
                        ranks = "matrix",
                        sens = "matrix",
                        spec = "matrix",
                        pAUC = "numeric",
                        AUC = "numeric",
                        factor = "factor",
                        cutpoints = "matrix",
                        caseNames = "character",
                        p = "numeric"),
         validity=function(object){
           if(any(dim(object@sens) != dim(object@spec)))
             return("\n'sens' and 'spec' must be matrices with equal dimensions")
           if(length(object@pAUC) != nrow(object@sens))
             return("\n'pAUC' must be numeric of length equal to nrow(sens)")
           if(length(object@factor)!=ncol(object@data) ||
              length(levels(object@factor))!=2)
             return("'factor' must be factor object with two levels and length = ncol(data)")
           if(length(object@pAUC) != length(object@AUC))
             return("'pAUC' and 'AUC' must be numeric vectors of equal length")
           if(nrow(object@cutpoints) != length(object@pAUC))
             return("'cutpoints' must be matrix with nrow=length(pAUC)")
           if(length(object@caseNames)!=2)
             return("'caseNames' must be character vector of length 2")
           return(TRUE)
         }
         )
## ==========================================================================
