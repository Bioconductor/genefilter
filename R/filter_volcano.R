filter_volcano <- function(
                           d, p, S,
                           n1, n2,
                           alpha, S_cutoff,
                           cex = .5, pch = 19,
                           xlab = expression( paste( log[2], " fold change" ) ),
                           ylab = expression( paste( "-", log[10], " p" ) ),
                           cols = c( "grey80", "grey50", "black" ),
                           ltys = c( 1, 3 ),
                           use_legend = TRUE,
                           ...
                           )
{

  f <- S < S_cutoff
  
  col <- rep( cols[1], length(d) )
  col[ !f & p >= alpha ] <- cols[2]
  col[ !f & p < alpha ] <- cols[3]  

  plot(
       d,
       -log10( p ),
       cex = cex,
       pch = pch,
       xlab = xlab,
       ylab = ylab,
       col = col,
       ...
       )

  k_grid <- seq( 0, max( -log10( p ) ), length = 100 )
  p_grid <- 10^( -k_grid )
  
  lines( kappa_p( p_grid, n1, n2 ) * S_cutoff, k_grid, lty = ltys[1] )
  lines( -1 * kappa_p( p_grid, n1, n2 ) * S_cutoff, k_grid, lty = ltys[1] )

  segments(
           c( par("usr")[1], kappa_p( alpha, n1, n2 ) * S_cutoff ),
           -log10( alpha ),
           c( -kappa_p( alpha, n1, n2 ) * S_cutoff, par("usr")[2] ),
           -log10( alpha ),
           lty = ltys[2]
           )

  if ( use_legend )
    legend(
           "topleft",
           c( "Filtered", "Insig.", "Sig." ),
           pch = pch,
           col = cols,
           inset = .025,
           bg = "white"
           )

}

