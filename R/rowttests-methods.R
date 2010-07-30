##--------------------------------------------------------------------------------------##
## This file contains methods definitions for rowttests, colttest, rowFtests, colFtests ##
##--------------------------------------------------------------------------------------##


##-----------------------------------------------------------------------
## The core function for row- and column-wise t-tests - it uses C code
##------------------------------------------------------------------------
rowcoltt =  function(x, fac, tstatOnly, which) {
  if (!missing(tstatOnly) && (!is.logical(tstatOnly) || is.na(tstatOnly)))
    stop(sQuote("tstatOnly"), " must be TRUE or FALSE.")
  
  f = checkfac(fac)
  if ((f$nrgrp > 2) || (f$nrgrp <= 0))
    stop("Number of groups is ", f$nrgrp, ", but must be >0 and <=2 for 'rowttests'.")

  cc = .Call("rowcolttests", x, f$fac, f$nrgrp, which-1L, PACKAGE="genefilter")
    
  res = data.frame(statistic = cc$statistic,
                   dm        = cc$dm,
                   row.names = dimnames(x)[[which]])

  if (!tstatOnly)
    res = cbind(res, p.value = 2*pt(abs(res$statistic), cc$df, lower.tail=FALSE))

  attr(res, "df") = cc$df    
  return(res)
}

##------------------------------------------------------------
## The core function for F-tests - it uses R matrix algebra 
##------------------------------------------------------------
rowcolFt =  function(x, fac, var.equal, which) {
  
  if(!(which %in% c(1L, 2L)))
    stop(sQuote("which"), " must be 1L or 2L.")
  
  if(which==2L)
    x = t(x)

  sqr = function(x) x*x
  
  stopifnot(length(fac)==ncol(x), is.factor(fac), is.matrix(x))
  x   <- x[,!is.na(fac), drop=FALSE]
  fac <- fac[!is.na(fac)]

  ## Number of levels (groups)
  k <- nlevels(fac)

  ## xm: a nrow(x) x nlevels(fac) matrix with the means of each factor
  ## level
  xm <- matrix(sapply(levels(fac), function(fl) rowMeans(x[,which(fac==fl), drop=FALSE])),nrow=nrow(x))

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
  
  res = data.frame(statistic = fstat,
                   p.value   = pf(fstat, dff, dfr, lower.tail=FALSE),
                   row.names = rownames(x))

  attr(res, "df") = c(dff=dff, dfr=dfr)
  return(res)
}

## ==========================================================================
## rowttests and colttests methods for 'matrix'
## ==========================================================================
setMethod("rowttests", signature(x="matrix", fac="factor"),
          function(x, fac, tstatOnly=FALSE)
          rowcoltt(x, fac, tstatOnly, 1L))

setMethod("rowttests", signature(x="matrix", fac="missing"),
          function(x, fac, tstatOnly=FALSE) 
          rowcoltt(x, factor(integer(ncol(x))), tstatOnly, 1L))

setMethod("colttests", signature(x="matrix", fac="factor"),
          function(x, fac, tstatOnly=FALSE)
          rowcoltt(x, fac, tstatOnly, 2L))

setMethod("colttests", signature(x="matrix", fac="missing"),
          function(x, fac, tstatOnly=FALSE) 
          rowcoltt(x, factor(integer(ncol(x))), tstatOnly, 2L))


## ==========================================================================
## rowFtests and colFtests methods for 'matrix'
## ==========================================================================
setMethod("rowFtests", signature(x="matrix", fac="factor"),
          function(x, fac, var.equal=TRUE)
          rowcolFt(x, fac, var.equal, 1L))

setMethod("colFtests", signature(x="matrix", fac="factor"),
          function(x, fac, var.equal=TRUE)
          rowcolFt(x, fac, var.equal, 2L))


## ===========================================================================
## Methods for 'ExpressionSet': only for rowttests and rowFtests
## -==========================================================================
setMethod("rowttests", signature(x="ExpressionSet", fac="factor"),
  function(x, fac, tstatOnly=FALSE)
    rowcoltt(exprs(x), fac, tstatOnly=tstatOnly, 1L))

setMethod("rowttests", signature(x="ExpressionSet", fac="missing"),
  function(x, fac, tstatOnly=FALSE) {
    x = exprs(x)
    fac = integer(ncol(x))
    rowcoltt(x, fac, tstatOnly, 1L)
  })

setMethod("rowttests", signature(x="ExpressionSet", fac="character"),
  function(x, fac, tstatOnly=FALSE) {
    if (length(fac) != 1)
      stop("fac must be length 1 character or a factor")
    fac = factor(pData(x)[[fac]])
    rowcoltt(exprs(x), fac, tstatOnly, 1L)
  })

setMethod("rowFtests", signature(x="ExpressionSet", fac="factor"),
  function(x, fac, var.equal=TRUE)
    rowcolFt(exprs(x), fac, var.equal, 1L))

setMethod("rowFtests", signature(x="ExpressionSet", fac="character"),
 function(x, fac, var.equal=TRUE) {
   fac = factor(as.integer(factor(pData(x)[[fac]]))-1L)
   rowcolFt(exprs(x), fac, var.equal, 1L)
 })



## ------------------------------------------------------------
## convert fac from factor or numeric to integer and then
## make sure it is an integer 
## ------------------------------------------------------------
checkfac = function(fac) {

  if(is.numeric(fac)) {
    nrgrp = as.integer(max(fac, na.rm=TRUE)+1)
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
  
  if(any(fac<0, na.rm=TRUE))
    stop("'fac' must not be negative.")
    
  return(list(fac=fac, nrgrp=nrgrp))
}




