shorth <- function(x, na.rm=FALSE) {
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
    if (max(q)-min(q) <= 0.05*length(x)) {
      q <- mean(q)
    } else {
      stop(paste("Encountered a tie in the calculation of q: ", q, "\n",
                 sx[q], "\n", sx[q+width], "\n"))
    }
    rv <- mean(sx[q:(q+width-1)])
  }
  return(rv)
}

