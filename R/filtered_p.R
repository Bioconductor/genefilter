filtered_p <- function( filter, test, theta, data, method = "none" ) {

  if ( is.function( filter ) )
    U1 <- filter( data )
  else
    U1 <- filter

  cutoffs <- quantile( U1, theta )

  result <- matrix( NA_real_, length( U1 ), length( cutoffs ) )
  colnames( result ) <- names( cutoffs )
  
  for ( i in 1:length( cutoffs ) ) {    
    use <- U1 >= cutoffs[i]
    if( any( use ) ) {
      if( is.function( test ) )
        U2 <- test( data[use,] )
      else
        U2 <- test[use]
      result[use,i] <- p.adjust( U2, method )
    }
  }

  return( result )
  
}




filtered_R <- function( alpha, filter, test, theta, data, method = "none" ) {

  p <- filtered_p( filter, test, theta, data, method )

  return( apply( p, 2, function(x) sum( x < alpha, na.rm = TRUE ) ) )

}
