
# genefinder.R
#
# genefinder functions.


scaleVector <- function(v) {
#
# Scales the elements of the vector v,
# and returns the result as a vector.
#
    mm <- range(v)
    result <- (v - mm[1]) / (mm[2] - mm[1])
    result
}



genescale <- function (m, axis=2) {
#
# scales a matrix using the scaleVector function.
#
   if( is.matrix(m) || is.data.frame(m) ) {
       rval <- apply (m, axis, scaleVector)
       if( axis==1 ) return(t(rval))
       return(rval)
   }
   else
	scaleVector(m)
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


    ninterest <- length(ilist)
    if( is.character(ilist) ) {
        iRows <- match(ilist, row.names(X))
        names(iRows) <- ilist
    }
    rval <- vector("list", length=ninterest)
    for (i in 1:ninterest) {
        ## xvals <- as.numeric( X[ilist[i], ] )
        ## rval[[i]] <- apply (X, 1, cor, xvals)

        N <- nrow(X)
        rval [[i]] <- .C("mm_distance",
	                 X = as.double(X),
                         nr= N,
	                 nc= ncol(X),
	                 d = double(N),
	                 iRow  = as.integer(ilist[i]),
	                 method= as.integer(method),
	                 DUP = FALSE, NAOK=TRUE, PACKAGE="genefilter")$d
    }
    return(rval)
}






