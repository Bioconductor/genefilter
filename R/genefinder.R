
# genefinder.R
#
# genefinder functions.


scaleVector <- function(v, na.rm=TRUE) {
#
# Scales the elements of the vector v,
# and returns the result as a vector.
#
    mm <- range(v, na.rm=na.rm)
    result <- (v - mm[1]) / (mm[2] - mm[1])
    result
}



genescale <- function (m, axis=2, na.rm=TRUE) {
#
# scales a matrix using the scaleVector function.
#
   if( is.matrix(m) || is.data.frame(m) ) {
       rval <- apply (m, axis, scaleVector, na.rm=na.rm)
       if( axis==1 ) return(t(rval))
       return(rval)
   }
   else
	scaleVector(m, na.rm=na.rm)
}



genefinder <- function (X, ilist, scale="none", method="euclidean") {
#
#
    X <- as.matrix(X)
    METHODS<- c("euclidean","maximum","manhattan","canberra","binary")
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
    ninterest <- length(ilist)
    if( is.character(ilist) ) {
        iRows <- match(ilist, row.names(X))
        names(iRows) <- ilist
    }
    else if ( is.numeric(ilist) )
        iRows <- ilist
    else
        stop("invalid genes selected")
    rval <- vector("list", length=ninterest)
    N <- nrow(X)
    for (i in 1:ninterest) {
        rval [[i]] <- .C("mm_distance",
	                 X = as.double(X),
                         nr= N,
	                 nc= ncol(X),
	                 d = double(N),
	                 iRow  = as.integer(iRows[i]),
	                 method= as.integer(method),
	                 DUP = FALSE, NAOK=TRUE, PACKAGE="genefilter")$d
    }
    return(rval)
}






