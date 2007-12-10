
##FIXME: this could replace the code further below at some point,
## but only when it has the var.equal option
##--------------------------------------------------
## fastT
##--------------------------------------------------
#fastT = function(x, ig1, ig2, var.equal=TRUE) {
#  fac      = rep(NA, ncol(x))
#  fac[ig1] = 0
#  fac[ig2] = 1
#  .Call("rowcolttests", x, as.integer(fac), as.integer(2),
#               as.integer(0), PACKAGE="genefilter")
#}



fastT = function(x, ig1, ig2, var.equal=TRUE) {
    ng1=length(ig1)
    ng2 = length(ig2)
    if( ncol(x) != ng1+ng2)
        stop("wrong sets of columns")

    outd = x[,c(ig1, ig2),drop=FALSE]
    nr = nrow(outd)
    z = rep(0, nr)
    dm = rep(0, nr)
    Z = .Fortran("fastt", d=as.single(outd), as.integer(nr),
           as.integer(ng1+ng2), as.integer(ng1), z = as.single(z),
         dm = as.single(dm), var.equal=as.integer(var.equal),
         ratio = as.integer(as.integer(0)), PACKAGE="genefilter")
    return(list(z = Z$z, dm=Z$dm, var.equal=Z$var.equal))
}
