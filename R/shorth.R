shorth <- function(x, na.rm=FALSE, tie.action="mean", tie.limit=0.05) {
  stopifnot(is.numeric(x))
  if (na.rm) {
    x <- x[is.finite(x)]
  } else {
    if(any(!is.finite(x)))
      stop("'x' contains NA or NaN, and 'na.rm' is FALSE.")
  }
    
  if(length(x)==0L) {
    
    NA_real_
      
    } else {
      
      sx    <- sort(x)
      width <- round(0.5*length(x))
      diffs <- sx[(width+1):length(x)] - sx[seq_len(length(x)-width)]
      
      ## cannot use which.min since we want all minimising points not just one:
      q  <- which(diffs==min(diffs))
      
      if(length(q)>1) {
        ## deal with ties:
        maxq = max(q)
        minq = min(q)
        ## take the action specified in "tie.action"
        q <- switch(tie.action,
                    mean = {
                      if (maxq-minq <= tie.limit * length(x)) {
                        mean(q)
                      } else {
                        stop(paste("Encountered tie(s), and the difference between minimal and maximal value is larger than 'length(x)*tie.limit'.",
                                   "This could mean that the distribution does not have a single well-defined mode.",
                                   paste("q=", minq, "...", maxq, ",  values=", signif(sx[minq],4), "...", signif(sx[minq+width],4), sep=""), sep="\n"))
                      }},
                    max  = maxq, ## largest midpoint (maxq)
                    min  = minq, ## smallest midpoint (minq)
                    stop(sprintf("Invalid value '%s' for argument 'tie.action'", tie.action))
                    )
      } ## if
      
      mean(sx[q:(q+width-1)])
      
    } ## if
  
}

