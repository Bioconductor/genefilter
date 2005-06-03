rowpAUCs <- function(x, fac, cutpts, p=0.1){

  fac = as.integer(fac)
  if (!all(sort(unique(fac)) == c(0, 1)))
        stop("'fac' variable must take values 0 or 1")

  ## TO DO: also cater for exprSets, as in rowttests
  if(!is.matrix(x))
    stop("'x' must be a matrix.")

  if(!is.numeric(p) || length(p)>1)
    stop("'p' must be numeric of length 1")
  
  # cutpoints
  if (missing(cutpts)) {
    doCuts = function(w){
      w = as.numeric(w)
      sw = sort(w)
      delta = min(diff(sw))/2
      cutpts = c(sw - delta, sw[length(sw)] + delta)
    }
    cutpts = t(apply(x, 1, doCuts))
  } else {
    if(length(cutpts)<1)
      stop("invalid number of cutpoints")
    if(!is.matrix(cutpts))
      cutpts <- matrix(cutpts, ncol=length(cutpts), nrow=nrow(x), byrow=TRUE)
    else if(ncol(cutpts) != ncol(x))
      stop("invalid cutpoints matrix") 
  }

  #memory allocation 
  nd = as.integer(dim(x))
  nc = as.integer(ncol(cutpts))
  spec = sens = matrix(0, nrow=nd[1], ncol=nc)
  pAUC = numeric(nd[1])
  return(invisible(.C("pAUC", data=x, nd=nd, cutpts=cutpts, nc=nc,
            truth=fac, spec=spec, sens=sens, pAUC=pAUC,
            p=p))[6:8])
}
