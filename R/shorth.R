shorth <- function(x, na.rm=FALSE, tieLimit=0.05) {
  stopifnot(is.numeric(x))
  if (na.rm)
    x <- x[!is.na(x)]
  rv <- NA
  if(length(x)>=1) {
    sx    <- sort(x)
    width <- round(0.5*length(x))
    diffs <- sx[(width+1):length(x)] - sx[1:(length(x)-width)]
    q     <- which(diffs==min(diffs))

    ## deal with ties: if they lie within 5%, simply take the average
    ## otherwise, generate an error
    maxq = max(q)
    minq = min(q)
    if (maxq-minq <= tieLimit * length(x)) {
      q <- mean(q)
    } else {
      stop(paste("Encountered a tie, this could mean that the distribution does not have a well-defined peak.\nq=",
                 minq, "...", maxq, "\nvalues: ", signif(sx[minq],4), "...", signif(sx[minq+width],4), "\n", sep=""))
    }
    rv <- mean(sx[q:(q+width-1)])
  }
  return(rv)
}

