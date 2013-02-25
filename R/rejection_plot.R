rejection_plot <- function(p,
                           col, lty = 1, lwd = 1,
                           xlab = "p cutoff", ylab = "number of rejections",
                           xlim = c( 0, 1 ), ylim,
                           legend = names(p),
                           at = c( "all", "sample" ),
                           n_at = 100,
                           probability = FALSE,
                           ...
                           )
{

  if ( is.matrix( p ) ) {
    legend <- colnames( p )
    p <- lapply( 1:ncol(p), function(i) p[,i] )
  }

  if ( missing( col ) )
    col <- rainbow( length( p ), v = .7 )

  col <- rep( col, length.out = length( p ) )
  lty <- rep( lty, length.out = length( p ) )
  lwd <- rep( lwd, length.out = length( p ) )
    
  if ( missing( ylim ) )
    ylim <- c( 0, ifelse( probability, 1, max( sapply( p, length ) ) ) )

  at <- match.arg( at )
  
  steps <- lapply(
                  p,
                  function(x) {
                    x <- na.omit(x)
                    stepfun(
                            sort( x ),
                            ( 0:length(x) ) / ifelse( probability, length(x), 1 )
                            )
                  }
                  )

  plot(
       0,
       type = "n",
       xaxs = "i", yaxs = "i",
       xlim = xlim, ylim = ylim,
       xlab = xlab, ylab = ylab,
       ...
       )

  if ( at == "all" ) {
    for ( i in 1:length( steps ) )
      lines(
            steps[[i]],
            xlim = xlim,
            col = col[i], lty = lty[i], lwd = lwd[i],
            do.points = FALSE
            )
  }
  
  else {
    x <- seq( xlim[1], xlim[2], length = n_at )
    for ( i in 1:length( steps ) )
      lines(
            x, steps[[i]](x),
            col = col[i], lty = lty[i], lwd = lwd[i]
            )
  }

  if ( !is.null( legend ) )
    legend(
           "topleft", 
           legend,
           col = col, lty = lty, lwd = lwd,           
           inset = .05
           )
  
  invisible( steps )

}
