# genefinder.R
#
# genefinder functions.


genescale <- function (m, axis=2, method=c("Z", "R"), na.rm=TRUE) {
    ##scale by the range
    RscaleVector <- function(v, na.rm=TRUE) {
        mm <- range(v, na.rm=na.rm)
        (v - mm[1]) / (mm[2] - mm[1])
    }
    ##scale using Zscore
    ZscaleVector <- function(v, na.rm=TRUE)
        (v - mean(v, na.rm=na.rm))/sd(v, na.rm=na.rm)
#
# scales a matrix using the scaleVector function.
#
    which <- match.arg(method)
    method <- switch(which,
                     Z = ZscaleVector,
                     R = RscaleVector)
    if( is.matrix(m) || is.data.frame(m) ) {
        rval <- apply (m, axis, method, na.rm=na.rm)
        if( axis==1 ) return(t(rval))
        return(rval)
    }
    else
	method(m, na.rm=na.rm)
}

.initFinder <- function(where) {

setGeneric("genefinder", function(X, ilist, numResults=25, scale="none",
    method="euclidean" )
    standardGeneric("genefinder"), where=where)

setMethod("genefinder", c("exprSet", "vector", "ANY", "ANY", "ANY"),
          function(X, ilist, numResults, scale,
                        method) {
          ans <- genefinder(exprs(X), ilist, numResults, scale,
                        method=method)
          ans$names <- geneNames(X)[ans$indices]
          ans
      }, where=where)

setMethod("genefinder", c("matrix", "vector", "ANY", "ANY", "ANY"),
         function (X, ilist, numResults, scale,
                        method) {
#
#
    print(method)
    X <- as.matrix(X)
    METHODS<-c("euclidean", "maximum", "manhattan", "canberra",
                        "correlation", "binary")
    method<-pmatch(method,METHODS)
    if (is.na(method))
        stop ("The distance method is invalid.")

    # perform scaling if requested.
    #
    if (scale == "none") {
        # no scaling
    } else
    if (scale == "range") {
       # scale input matrix using 'genescale'
        X <- genescale(X)
    } else
    if (scale == "zscore") {
        # scale using R's scale function
        X <- scale(X)
    } else
    stop ("The scale method is invalid.")

    if( !is.vector(ilist) )
        stop("the genes to be compared to must be in a vector")

    ninterest <- length(ilist);

    if( is.character(ilist) ) {
        iRows <- match(ilist, row.names(X))
        names(iRows) <- ilist
    }
    else if ( is.numeric(ilist) )
        iRows <- ilist
    else
        stop("invalid genes selected")

    N <- nrow(X)

    Genes <- array(dim=c(ninterest, numResults))
    Dists <- array(dim=c(ninterest, numResults))
    extCall <- .C("mm_distance",
                  X = as.double(X),
                  nr= N,
                  nc= ncol(X),
                  g = as.integer(Genes),
                  d = as.double(Dists),
                  iRow  = as.integer(iRows),
                  nInterest = as.integer(ninterest),
                  nResults = as.integer(numResults),
                  method= as.integer(method),
                  DUP = FALSE, NAOK=TRUE, PACKAGE="genefilter")

    Genes <- extCall$g
    Dists <- extCall$d
    return(list(indices=Genes, dists=Dists))
}, where=where)

}






