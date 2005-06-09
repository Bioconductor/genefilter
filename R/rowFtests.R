##------------------------------------------------------------
## x: a numeric matrix with n rows ("genes") and d columns ("conditions")
## fac: a factor
##------------------------------------------------------------
rowFtests = function(x, fac,var.equal=TRUE) {
  if(is(x, "exprSet")) {
    if(is.character(fac))
      fac = pData(x)[[fac]]
    x   = exprs(x)
  }
   
  sqr = function(x) x*x
  
  stopifnot(length(fac)==ncol(x), is.factor(fac), is.matrix(x))
  x   <- x[,!is.na(fac), drop=FALSE]
  fac <- fac[!is.na(fac)]

  ## Number of levels (groups)
  k      <- nlevels(fac)

  ## xm: a nrow(x) x nlevels(fac) matrix with the means of each factor
  ## level
  xm <- sapply(levels(fac), function(fl) rowMeans(x[,which(fac==fl), drop=FALSE]))

  ## x1: a matrix of same size as x with group means
  x1 <- xm[,fac, drop=FALSE]

  ## degree of freedom 1
  dff    <- k - 1

  if(var.equal){
    ## x0: a matrix of same size as x with overall means
    x0 <- matrix(rowMeans(x), ncol=ncol(x), nrow=nrow(x))
  
    ## degree of freedom 2
    dfr    <- ncol(x) - dff - 1

    ## mean sum of squares
    mssf   <- rowSums(sqr(x1 - x0)) / dff
    mssr   <- rowSums(sqr( x - x1)) / dfr

    ## F statistic
    fstat  <- mssf/mssr

  } else{

    ## a nrow(x) x nlevels(fac) matrix with the group size  of each factor
    ## level
    ni <- t(matrix(tapply(fac,fac,length),ncol=nrow(x),nrow=k))

    ## wi: a nrow(x) x nlevels(fac) matrix with the variance * group size of each factor
    ## level
    sss <- sqr(x-x1)
    x5 <- sapply(levels(fac), function(fl) rowSums(sss[,which(fac==fl), drop=FALSE]))
    wi <- ni*(ni-1) /x5

    ## u : Sum of wi
    u  <- rowSums(wi)

    ## F statistic
    MR <- rowSums(sqr((1 - wi/u)) * 1/(ni-1))*1/(sqr(k)-1)
    fsno <- 1/dff * rowSums(sqr(xm - rowSums(wi*xm)/u) * wi)
    fsdeno <- 1+ 2* (k-2)*MR
    fstat <- fsno/fsdeno

    ## degree of freedom 2: Vector with length nrow(x)
    dfr <- 1/(3 * MR)
  
  }
  
  ## p-value for F-test
  pval   <- pf(fstat, dff, dfr, lower.tail = FALSE)
  return(list(statistic=fstat, p.value=pval, dof=c(dff,dfr)))
}

##--------------------------------------------------
## colFtests
##--------------------------------------------------
colFtests = function(x, fac,var.equal=TRUE) 
  rowFtests(t(x), fac, var.equal) 

##--------------------------------------------------
## rowttests
##--------------------------------------------------
rowttests = function(x, fac, tstatOnly=FALSE) {
  if(is(x, "exprSet")) {
    if(is.character(fac))
      fac = as.integer(factor(pData(x)[[fac]]))-1
    x   = exprs(x)
  }
  
  f   = checkfac(fac)
  if(f$nrgrp>2)
    stop("Number of groups must be <= 2 for 'rowttests'.")
  
  res = .Call("rowcolttests", x, f$fac, f$nrgrp,
               as.integer(0), PACKAGE="genefilter")
  if(!tstatOnly)
    res$p.value = 2*pt(abs(res$statistic), res$df, lower.tail = FALSE)
  return(res)
}

##--------------------------------------------------
## colttests
##--------------------------------------------------
colttests = function(x, fac, tstatOnly=FALSE) {
  f   = checkfac(fac)
  res = .Call("rowcolttests", x, f$fac, f$nrgrp,
               as.integer(1), PACKAGE="genefilter")
  if(!tstatOnly)
    res$p.value = 2*pt(abs(res$statistic), res$df, lower.tail = FALSE)
  return(res)
}

##FIXME: this should replace the code below at some point,
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
#
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

## ------------------------------------------------------------
## convert fac from factor or numeric to integer and then
## make sure it is an integer 
## ------------------------------------------------------------
checkfac = function(fac) {
  if(missing(fac)) {
    fac   = integer(ncol(x))
  }
  if(is.numeric(fac)) {
    nrgrp = as.integer(max(fac)+1)
    fac   = as.integer(fac)
  }
  ## this must precede the factor test
  if(is.character(fac))
    fac = factor(fac)

  if (is.factor(fac)) {
    nrgrp = nlevels(fac)
    fac   = as.integer(as.integer(fac)-1)
  } 
  if(!is.integer(fac))
    stop("'fac' must be factor, character, numeric, or integer.")
  if(any(fac<0))
    stop("'fac' must not be negative.")
    
  return(list(fac=fac, nrgrp=nrgrp))
}
