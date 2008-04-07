## Generic functions for package genefilter

setGeneric("rowFtests", function(x, fac, var.equal=TRUE)
           standardGeneric("rowFtests"))
setGeneric("colFtests", function(x, fac, var.equal=TRUE)
           standardGeneric("colFtests"))
setGeneric("rowttests", function(x, fac, tstatOnly=FALSE)
           standardGeneric("rowttests"))
setGeneric("colttests", function(x, fac, tstatOnly=FALSE)
           standardGeneric("colttests"))


setGeneric("genefinder", function(X, ilist, numResults=25, scale="none",
    weights, method="euclidean" )
    standardGeneric("genefinder"))


setGeneric("pAUC", function(object, p, flip=TRUE) standardGeneric("pAUC"))

setGeneric("AUC", function(object) standardGeneric("AUC"))

setGeneric("sens", function(object) standardGeneric("sens"))

setGeneric("spec", function(object) standardGeneric("spec"))

setGeneric("area", function(object, total=FALSE) standardGeneric("area"))

setGeneric("rowpAUCs", function(x, fac, p=0.1, flip=TRUE, caseNames=c("1", "2"))
           standardGeneric("rowpAUCs"))

setGeneric("nsFilter", signature="eset",
           function(eset,
                    require.entrez=TRUE,
                    require.GOBP=FALSE,
                    require.GOCC=FALSE,
                    require.GOMF=FALSE,
                    remove.dupEntrez=TRUE,
                    var.func=IQR, var.cutoff=0.5, var.filter=TRUE,
                   filterByQuantile=TRUE,
                    feature.exclude="^AFFX", ...)
           standardGeneric("nsFilter"))

