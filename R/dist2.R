dist2 = function(x, fun=function(a,b) median(abs(a-b), na.rm=TRUE), diagonal=0) {
  if(!(is.numeric(diagonal)&&(length(diagonal)==1L)))
    stop("'diagonal' must be a numeric scalar.")
  res = matrix(diagonal, ncol=ncol(x), nrow=ncol(x))
  colnames(res) = rownames(res) = colnames(x)
  if(ncol(x)>=2) {
    for(j in 2:ncol(x))
      for(i in 1:(j-1))
        res[i, j] = res[j, i] = fun(x[,i], x[,j])
  }
  return(res)
}
