##------------------------------------------------------------
## x: a numeric matrix with n rows ("genes") and d columns ("conditions")
## fac: a factor
##------------------------------------------------------------
rowFtests <- function(x, fac,var.equal=TRUE) {
   
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
## 
##--------------------------------------------------
rowttests = function(x, fac, tstatOnly=FALSE) {
  if(missing(fac)) {
    ## one group
    sqr   = function(x) x*x
    xavg  = rowMeans(x)
    denom = sqrt(rowSums(sqr(x-xavg)) / ((ncol(x)-1)*ncol(x)) )
    res   = list(statistic=xavg/denom, df=ncol(x)-1)
  } else {
    ## two groups
    fac = checkfac(fac)    
    res = .Call("rowcolttests", x, fac, as.integer(0), PACKAGE="genefilter")
  }
  if(!tstatOnly)
    res$p.value = 2*pt(abs(res$statistic), res$df, lower.tail = FALSE)
  return(res)
}

## ------------------------------------------------------------
## convert fac from factor or numeric to integer and then
## make sure it is an integer 
## ------------------------------------------------------------
checkfac = function(fac) {
  if(is.factor(fac))
    fac = as.integer(fac)-1;
  if(is.numeric(fac))
    fac = as.integer(fac);
  if(!is.integer(fac))
    stop("'fac' must be factor, numeric, or integer")
  return(fac)
}
