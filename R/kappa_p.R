kappa_p <- function( p, n1, n2 = n1 ) {
  n <- n1 + n2
  t <- qt( 1 - p/2, df = n - 2 )
  kappa_t( t, n1, n2 )
}

kappa_t <- function( t, n1, n2 = n1 ) {
  n <- n1 + n2
  sqrt( n * (n-1) * t^2 / ( n1 * n2 * ( n - 2 + t^2 ) ) )
}
