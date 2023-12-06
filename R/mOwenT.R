mOwenT <- function( h, r, ph, prh, atanExt = TRUE, xy = FALSE ) {
  # The function calculates the Owen''s function T( h, r ). It can handle the
  # special cases r = Inf and r = -Inf. In such cases, indefinite expressions
  # occur but are caught and modified so that the loop stops with s = 0, and
  # then the final result is formed. The calling function ensures that h and
  # r are of the same length and provides pnorm( h ) and pnorm( r * h ) in ph
  # and prh parameters respectively, where prh must be set to 0 if h is zero
  # and r is infinite. If xy = TRUE, the two calculations needed to calculate
  # Phi2xy( x, y ) are combined into one, since the calling function ensures
  # that h = c( x, y ), r = c( rx, ry ) etc.
  if ( xy ) {
    kx <- 1 : ( length( h ) / 2 )
    ky <- ( length( h ) / 2 + 1 ) : length( h )
  }
  u <- ph
  i <- abs( r ) > 1 # indicates where the transformations must be done
  # h and r are not transformed since it is easy to calculate p, q and a
  u[ i ] <- ( u / 2 + prh * ( 1 / 2 - u ) - ( 1 - sign( r ) ) / 4 )[ i ]
  # Parameter checking in the calling function ensures that h and r have the
  # same precBits if they are mpfr numbers. If Rmpfr is used, all variables
  # below are automatically given class mpfr due to calculation with h or r,
  # or with other already transformed variables.
  if ( isa( h, "mpfr" ) ) pi <- Rmpfr::Const( "pi", Rmpfr::getPrec( h ) )
  p <- pmin( r^2, 1 ) / ( 1 + r^2 )
  q <- ( 1 + r^2 ) * h^2 / 2 # q is invariant under parameter transformation
  a <- abs( r ) / ( 2 * pi * ( 1 + r^2 ) ) # a is also invariant
  b <- exp( -q )
  d <- b
  if ( atanExt ) s <- a * ( 1 - d ) else s <- a * d
  # n is number of iterations needed for the components
  if ( xy ) n <- rep( 0, length( s ) / 2 ) else n <- rep( 0, length( s ) )
  k <- 0
  repeat {
    a <- ( 2 * k + 2 ) * p * a / ( 2 * k + 3 )
    b <- q * b / ( k + 1 )
    d <- d + b
    if ( atanExt ) z <- s + a * ( 1 - d ) else z <- s + a * d
    z[ is.nan ( z ) ] <- 0 # needed because of special cases
    k <- k + 1
    if ( xy )
      n[ n == 0 & z[ kx ] + z[ ky ] <= s[ kx ] + s[ ky ] ] <- k
    else
      n[ n == 0 & z <= s ] <- k
    if ( all( n > 0 ) ) break
    s <- z
  }
  if ( atanExt ) {
    z <- atan( abs ( r ) ) / ( 2 * pi )
    z[ i ] <- 1 / 4 - z[ i ] # = atan( abs ( 1 / r[ i ] ) ) / ( 2 * pi )
    s <- z - s
  }
  s[ r < 0 ] <- -s[ r < 0 ] # equating the sign with the sign of r
  s[ i ] <- u[ i ] - s[ i ] # adjustment due to parameter transformation
  if ( xy ) s <- s[ kx ] + s[ ky ]
  attr( s, "nIter" ) <- n
  return( s )
}
