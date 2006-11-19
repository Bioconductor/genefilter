binaryRepresentation = function(n, ndigits) {

  if(any(n<0))
    stop("Please do not specify negative numbers in 'n'")
  
  ndig = max(floor(log2(n)))+as.integer(1)
  if(missing(ndigits)){
    ndigits=ndig
  } else {
    if(any(ndig>ndigits))
      stop(sprintf("Numbers in 'n' require at least %d digits, but you specified only ndigits=%d", ndig, ndigits)) 
  }
  
  res = matrix(NA, nrow=ndigits, ncol=length(n))
  expt = 2^(0:(ndigits-1))
  colnames(res) = paste(n)
  rownames(res) = paste(expt)
  for(j in ndigits:1) {
    k = (n %/% expt[j]) > 0
    res[j,] = k
    n = n - expt[j] * as.integer(k)
  }
  if(any(n!=0))
    stop("Internal error, please give note to function author.")
  
  return(res)
}

