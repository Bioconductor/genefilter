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
  return(list(statistic=fstat, pvalue=pval, dof=c(dff,dfr)))
}

rowttests <- function(x, fac) {
  stopifnot(is.matrix(x))
  sqr <- function(z) z*z

  if(missing(fac)) {
    ## one group
    xavg  <- rowMeans(x)
    dof   <- ncol(x)-1
    denom <- sqrt(rowSums(sqr(x-xavg)) / (dof*ncol(x)) )
    tstat <- xavg/denom
  } else {
    ## two groups
    if(length(fac)!=ncol(x) || !is.factor(fac))
      stop("'fac' must be a factor with length 'ncol(x)'.")
    levs <- levels(fac)
    if(!length(levs)==2)
      stop("'fac' must have exactly two levels.")
    x   <- x[,!is.na(fac), drop=FALSE]
    fac <- fac[!is.na(fac), drop=FALSE]
  
    g1 <- which(fac==levs[1])    
    g2 <- which(fac==levs[2])
    n1 <- length(g1)
    n2 <- length(g2)
    xavg1 <- rowMeans(x[,g1, drop=FALSE])
    xavg2 <- rowMeans(x[,g2, drop=FALSE])
    dof <- n1+n2-2
    denom <- sqrt( (rowSums(sqr(x[,g1, drop=FALSE] - xavg1)) +
                    rowSums(sqr(x[,g2, drop=FALSE] - xavg2))) * ((n1+n2)/(dof*n1*n2)) )
    tstat <- (xavg2 - xavg1) / denom
  }
  pval  <- 2*pt(abs(tstat), dof, lower.tail = FALSE)
  return(list(statistic=tstat, pvalue=pval, dof=dof))
}

