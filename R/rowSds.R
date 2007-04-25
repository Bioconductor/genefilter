rowVars = function(x, ...) {
  sqr     = function(x)  x*x
  n       = rowSums(!is.na(x))
  n[n<=1] = NA
  return(rowSums(sqr(x-rowMeans(x, ...)), ...)/(n-1))
}

rowSds = function(x, ...)
  sqrt(rowVars(x, ...))
