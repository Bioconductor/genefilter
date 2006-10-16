doCuts = function(x){
  sx = sort(as.numeric(x))
  delta = min(diff(sx))/2
  cutpts = c(sx - delta, sx[length(sx)] + delta)
}


rowpAUCs <- function(x, fac, cutpts, p=0.1){
    ## FIXME: this should become generic and method
  if(is(x, "exprSet") || is(x, "ExpressionSet")) {
     if(is.character(fac))
       fac = as.integer(factor(pData(x)[[fac]]))-1
     x   = exprs(x)
  }
  f = genefilter:::checkfac(fac)
  if(f$nrgrp != 2 || length(f$fac) != ncol(x) || length(unique(f$fac)) !=2 )
    stop("'fac' must be factor with 2 levels and length 'ncol(x)'")

  if(!is.matrix(x))
    stop("'x' must be a matrix.")

  if(!is.numeric(p) || length(p)>1)
    stop("'p' must be numeric of length 1")
  
  # cutpoints
  if (missing(cutpts)) {
    cutpts = t(apply(x, 1, doCuts))
  } else {
    if(length(cutpts)<2)
      stop("invalid number of cutpoints")
    if(!is.matrix(cutpts))
      cutpts <- matrix(cutpts, ncol=length(cutpts), nrow=nrow(x), byrow=TRUE)
    if(ncol(cutpts) != ncol(x))
      stop("invalid cutpoints matrix") 
  }

  #memory allocation 
  #nd = as.integer(dim(x))
  #nc = as.integer(ncol(cutpts))
  #spec = sens = matrix(0, nrow=nd[1], ncol=nc)
  #pAUC = numeric(nd[1])
  #return(invisible(.C("pAUC", data=x, nd=nd, cutpts=cutpts, nc=nc,
  #          truth=as.integer(f$fac), spec=spec, sens=sens, pAUC=pAUC,
  #          p=p))[6:8])
  res <- .Call("pAUC", x, cutpts, as.integer(f$fac), p, PACKAGE="genefilter")
  return(res)
}
