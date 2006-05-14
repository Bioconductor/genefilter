half.range.mode <- function( data, B, B.sample, beta = .5, diag = FALSE ) {

  if (missing( B ) ) {

    # Just one run on the full set...
    
    if ( is.unsorted( data ) ) data <- sort( data )
    
    .C(
       "half_range_mode",
       data = as.double( data ),
       n = as.integer( length( data ) ),
       beta = as.double( beta ),
       diag = as.integer( diag ),
       M = double(1),
       PACKAGE = "symp"
       )$M

  }

  else {

    # Bootstrapped

    if ( missing( B.sample ) )
      B.sample <- length( data )

    M <- sapply(
                1:B,
                function (x) {
                  d <- sort( sample( data, B.sample, replace = T ) )
                  .C(
                     "half_range_mode",
                     data = as.double( d ),
                     n = as.integer( B.sample ),
                     beta = as.double( beta ),
                     diag = as.integer( diag ),
                     M = double(1),
                     PACKAGE = "symp"
                     )$M
                }
                )

    mean( M )
    
  }

}

