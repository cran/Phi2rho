vOwenT <- function( h, r, ph, prh, transf = TRUE, xy = FALSE ) {
  # The function calculates the Owen''s function T( h, r ). It can handle the
  # special cases r = Inf and r = -Inf, and also r = 0 if transf = TRUE. In
  # such cases, indefinite expressions occur but are caught and modified so
  # that the loop stops with Q = 0, and then the final result is formed. The
  # calling function ensures that h and r are of the same length and provides
  # pnorm( h ) and pnorm( r * h ) in ph and prh parameters respectively, where
  # prh must be set to 0 if h is zero and r is infinite. If xy = TRUE, the two
  # calculations needed to calculate Phi2xy( x, y ) are combined into one,
  # since the calling function ensures that h = c( x, y ), r = c( rx, ry ) etc.
  if ( !transf ) stopifnot( all( r != 0 ) )
  if ( xy ) {
    kx <- 1 : ( length( h ) / 2 )
    ky <- ( length( h ) / 2 + 1 ) : length( h )
  }
  u <- ph
  v <- ph
  i <- transf & abs( r ) < 1 # indicates where the transformations must be done
  # attention: the transformation of u must be before the transformation of r,
  # because otherwise for r = Inf and r = -Inf we get the sign 0 instead of 1
  # and -1 repectively, because sign( 0 ) = 0 in R
  v[ i ] <- prh[ i ]
  u[ i ] <- ( u + v * ( 1 - u ) - ( 1 - sign( r ) ) / 4 )[ i ]
  h[ i ] <- r[ i ] * h[ i ]
  r[ i ] <- 1 / r[ i ]
  # Parameter checking in the calling function ensures that h and r have the
  # same precBits if they are mpfr numbers. If Rmpfr is used, all variables
  # below are automatically given class mpfr due to calculation with h or r,
  # or with other already transformed variables.
  if ( isa( h, "mpfr" ) ) pi <- Rmpfr::Const( "pi", Rmpfr::getPrec( h ) )
  p <- 1 / ( 1 + r^2 )
  B <- sqrt( p ) * exp( -h^2 / ( 2 * p ) ) / ( 2 * pi )
  w <- -abs( h ) / sqrt( p ) # can be 0 / 0 if h is zero and r is infinite
  w[ is.nan( w ) ] <- 0 # needed because of special cases if Rmpfr is used
  A <- -abs( h ) * pnorm( w ) / sqrt( 2 * pi ) + B
  Q <- A
  # n is number of iterations needed for the components
  if ( xy ) n <- rep( 0, length( Q ) / 2 ) else n <- rep( 0, length( Q ) )
  k <- 0
  repeat {
    k <- k + 1
    B <- ( 2 * k - 1 )^2 * p * B / ( 2 * k * ( 2 * k + 1 ) )
    A <- -A * h^2 * ( 2 * k - 1 ) / ( 2 * k * ( 2 * k + 1 ) ) + B
    z <- Q + A
    z[ is.nan( z ) ] <- 0 # needed because of special cases
    if ( xy ) # == instead of <= because the series may not be monotone
      n[ n == 0 & z[ kx ] + z[ ky ] == Q[ kx ] + Q[ ky ] ] <- k
    else
      n[ n == 0 & z == Q ] <- k
    if ( all( n > 0 ) ) break
    Q <- z
  }
  Q[ r > 0 ] <- ( pmin( v, 1 / 2 ) - Q )[ r > 0 ]
  Q[ r < 0 ] <- ( pmax( v - 1 / 2, 0 ) + Q )[ r < 0 ]
  Q[ i ] <- -Q[ i ] + u[ i ] # adjustment due to parameter transformation
  Q <- Q - ph / 2 # adjustment to get T( h, r ) from Phi2( h, 0, rho )
  Q[ i & r == Inf ] <- 0 # needed because of special cases
  if ( xy ) Q <- Q[ kx ] + Q[ ky ]
  attr( Q, "nIter" ) <- n
  return( Q )
}
