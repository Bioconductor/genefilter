# genefinder.R
#
# genefinder functions.


genescale <- function (m, axis=2, method=c("Z", "R"), na.rm=TRUE) {
    ##scale by the range
    RscaleVector <- function(v, na.rm) {
        mm <- range(v, na.rm=na.rm)
        (v - mm[1]) / (mm[2] - mm[1])
    }
    ##scale using Zscore
    ZscaleVector <- function(v, na.rm)
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


setMethod("genefinder", c("ExpressionSet", "vector", "ANY", "ANY", "ANY",
          "ANY"),
          function(X, ilist, numResults, scale, weights,
                   method) {
              gN <- featureNames(X)
              if (is.character(ilist))
                  ilist <- match(ilist,gN)
              ans <- genefinder(exprs(X), ilist, numResults, scale, weights,
                        method=method)
              names(ans) <- gN[ilist]
              ans
      })

setMethod("genefinder", c("matrix", "vector", "ANY", "ANY", "ANY", "ANY"),
         function (X, ilist, numResults, scale, weights,
                        method) {
    X <- as.matrix(X)
    METHOD <- c("euclidean", "maximum", "manhattan",
                        "canberra", "correlation", "binary")
    method<-pmatch(method, METHOD)
    if (is.na(method))
        stop ("The distance method is invalid.")

    SCALE <- c("none", "range", "zscore")
    scale <- SCALE[pmatch(scale, SCALE)]

    # perform scaling if requested.
    #
    X <- switch(scale,
                none=X,
                range=genescale(X),
                zscore=scale(X),
                stop("The scaling method is invalid")
                )
    N <- nrow(X)
    C <- ncol(X)

    if( !is.vector(ilist) )
        stop("the genes to be compared to must be in a vector")

    ninterest <- length(ilist);

    if( is.character(ilist) ) {
        iRows <- match(ilist, row.names(X))
        names(iRows) <- ilist
    }
    else if ( is.numeric(ilist) ) {
        iRows <- ilist
        names(iRows) <- row.names(X)[ilist]
    }
    else
        stop("invalid genes selected")

    if( any(is.na(iRows)) )
        stop("invalid genes selected")

    if (missing(weights))
        weights <- rep(1,C)
    else if (length(weights) != C)
        stop("Supplied weights do not match number of columns")

    ## Do a sanity check on the requested genes in ilist -> if the
    ## gene exceeds the # of rows in the matrix, can not be processed.
    if (max(iRows) > N)
        stop("Requested genes exceed the dimensions of the supplied matrix.")


    Genes <- array(as.integer(NA), dim=c(ninterest, numResults))
    Dists <- array(as.integer(NA), dim=c(ninterest, numResults))
    extCall <- .C("gf_distance",
                  X = as.double(X),
                  nr= as.integer(N),
                  nc= ncol(X),
                  g = as.integer(Genes),
                  d = as.double(Dists),
                  iRow  = as.integer(iRows),
                  nInterest = as.integer(ninterest),
                  nResults = as.integer(numResults),
                  method= as.integer(method),
                  weights = as.double(weights),
                  DUP = FALSE, NAOK=TRUE, PACKAGE="genefilter")

    Genes <- extCall$g+1
    Dists <- extCall$d
    Which <- vector()

    ## Get the number of genes/dists per selection.  There should
    ## always be a number of total genes such that they are a multiple
    ## of ninterest
    numPerList <- length(Genes) / ninterest

    Which <- rep(iRows, rep(numPerList, ninterest))

    byGene <- split(Genes, Which)
    names(byGene) <- rep("indices", length(byGene))
    byDists <- split(Dists, Which)
    names(byDists) <- rep("dists", length(byDists))
    ## Need a better way to stuff these together
    retList <- list()
    for (i in 1:ninterest) {
        retList[[i]] <- list(indices=byGene[[i]], dists=byDists[[i]])
    }

    return(retList)
})







