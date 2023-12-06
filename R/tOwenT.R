tOwenT <- function( h, r, ph, prh, transf = TRUE, xy = FALSE ) {
  # The function calculates the Owen''s function T( h, r ). It can handle the
  # special cases r = Inf and r = -Inf. In such cases, indefinite expressions
  # occur but are caught and modified so that the loop stops with s = 0, and
  # then the final result is formed. The calling function ensures that h and
  # r are of the same length and provides pnorm( h ) and pnorm( r * h ) in ph
  # and prh parameters respectively, where prh must be set to 0 if h is zero
  # and r is infinite. If xy = TRUE, the two calculations needed to calculate
  # Phi2xy( x, y ) are combined into one, since the calling function ensures
  # that h = c( x, y ), r = c( rx, ry ) etc. In the first published version on
  # CRAN, the calculation for transf = FALSE and infinite r was not possible.
  # For such cases, exceptionally, we calculate as if transf = TRUE, which is
  # also taken into account by the calling functions when preparing the ph and
  # prh parameters.
  if ( xy ) {
    kx <- 1 : ( length( h ) / 2 )
    ky <- ( length( h ) / 2 + 1 ) : length( h )
  }
  u <- ph 
  # i indicates where the transformations must be done
  i <- ( transf & abs( r ) > 1 ) | is.infinite( r )
  # attention: the transformation of u must be before the transformation of r,
  # because otherwise for r = Inf and r = -Inf we get the sign 0 instead of 1
  # and -1 repectively, because sign( 0 ) = 0 in R
  u[ i ] <- ( u / 2 + prh * ( 1 / 2 - u ) - ( 1 - sign( r ) ) / 4 )[ i ]
  h[ i ] <- r[ i ] * h[ i ]
  r[ i ] <- 1 / r[ i ]
  # Parameter checking in the calling function ensures that h and r have the
  # same precBits if they are mpfr numbers. If Rmpfr is used, all variables
  # below are automatically given class mpfr due to calculation with h or r,
  # or with other already transformed variables.
  if ( isa( h, "mpfr" ) ) pi <- Rmpfr::Const( "pi", Rmpfr::getPrec( h ) )
  He0 <- rep( 1, length( h ) ) # = He_0(h)
  He1 <- h                     # = He_1(h)
  p <- r^2 / ( 1 + r^2 )
  q <- sqrt( p ) * dnorm( h ) / sqrt( 2 * pi )
  s <- q
  z <- rep( 0, length( s ) )
  # n is number of iterations needed for the components
  if ( xy ) n <- rep( 0, length( s ) / 2 ) else n <- z
  k <- 0
  repeat {
    k <- k + 1
    He0 <- h * He1 - ( 2 * k - 1 ) * He0 # = He_{2k} / ( 2^k * k! )
    He1 <- h * He0 - 2 * k * He1 # = He_{2k+1} / ( 2^k * k! )
    He0 <- He0 / ( 2 * k )
    He1 <- He1 / ( 2 * k )
    q <- -q * p
    v <- s + He0 * q / ( 2 * k + 1 )
    v[ is.nan( v ) ] <- 0 # needed because of special cases
    # if h is a zero of the 2k-th Hermite polynomial, the new approximation is
    # equal to the last one, but the recursion must not be terminated, so the
    # recursion is terminated when the new approximation is equal to the last
    # one and penultimate one
    if ( xy ) { # == instead of <= because the series may not be monotone
      t1 <- v[ kx ] + v[ ky ] == s[ kx ] + s[ ky ]
      t2 <- s[ kx ] + s[ ky ] == z[ kx ] + z[ ky ]
      n[ n == 0 & t1 & t2 ] <- k
    } else n[ n == 0 & v == s & s == z ] <- k 
    if ( all( n > 0 ) ) break
    z <- s
    s <- v
  }
  s[ r < 0 ] <- -s[ r < 0 ] # equating the sign with the sign of r
  s[ i ] <- u[ i ] - s[ i ] # adjustment due to parameter transformation
  if ( xy ) s <- s[ kx ] + s[ ky ]
  attr( s, "nIter" ) <- n
  return( s )
}
