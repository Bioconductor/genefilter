rowSds <- function(x, ...) {
  sqr     = function(x)  x*x
  n       = rowSums(!is.na(x))
  n[n<=1] = NA
  return(sqrt(rowSums(sqr(x-rowMeans(x, ...)), ...)/(n-1)))
}
