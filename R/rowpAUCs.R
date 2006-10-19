doCuts = function(x, unique=FALSE){
  if(unique)
    sx = sort(unique(as.numeric(x)))
  else
    sx = sort(as.numeric(x))
  delta = min(diff(sx))/2
  cutpts = c(sx - delta, sx[length(sx)] + delta)
  return(cutpts)
}



rowpAUCs <- function(x, fac, cutpts, p=0.1, common=FALSE){
  if(is(x, "exprSet") | is(x, "eSet")) {
     if(is.character(fac))
       fac = as.integer(factor(pData(x)[[fac]]))-1
     x   = exprs(x)
  }
  f = genefilter:::checkfac(fac)
  
  if(!is.matrix(x))
    stop("'x' must be a matrix.")
  
  if(f$nrgrp != 2 || length(f$fac) != ncol(x) || length(unique(f$fac)) !=2 )
    stop("'fac' must be factor with 2 levels and length 'ncol(x)'")


  if(!is.numeric(p) || length(p)>1)
    stop("'p' must be numeric of length 1")
  
  # cutpoints
  if (missing(cutpts)) {
    if(common){
      cp = unique(genefilter:::doCuts(x, TRUE))
    cutpts <- matrix(cp, ncol=length(cp), nrow=nrow(x), byrow=TRUE)
    }else
    cutpts <- t(apply(x, 1, genefilter:::doCuts)) 
  } else {
    if(length(cutpts)<2)
      stop("invalid number of cutpoints")
    if(!is.matrix(cutpts))
      cutpts <- matrix(cutpts, ncol=length(cutpts), nrow=nrow(x), byrow=TRUE)
    if(ncol(cutpts) != ncol(x))
      stop("invalid cutpoints matrix") 
  }

  
  res <- .Call("pAUC", x, cutpts, as.integer(f$fac), p, PACKAGE="genefilter")
  return(res)
}
