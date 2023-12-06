OwenT <- function( h, a, opt = TRUE, fun = c( "mOwenT", "tOwenT", "vOwenT" ) ) {
  # The OwenT function is a wrapper for the mOwenT, tOwenT and vOwenT functions,
  # which compute the Owen's function T( h, a ), but use r instead of the original
  # Owen's parameter a, which is preserved here.
  chkArgs1( h = h, a = a, opt = opt )
  if ( length( h ) < length( a ) ) h <- rep( h, length( a ) )
  if ( length( a ) < length( h ) ) a <- rep( a, length( h ) )
  z <- h # dummy - forced to have a mpfr class if Rmpfr is used
  n <- rep( 0, length( z ) ) # frame for the number of iterations
  z[ a == 0 ] <- 0 # needed if fun = "vOwenT" and opt = FALSE
  z[ is.na( a ) ] <- NA
  i <- a != 0 & !is.na( a )
  if ( any( i ) ) {
    h   <- h[ i ]
    a   <- a[ i ]
    j   <- !is.nan( a * h )
    fun <- match.arg( fun )   
    if ( fun == "mOwenT" ) j <- j & abs( a ) > 1
    if ( fun == "tOwenT" ) j <- j & ( opt & abs( a ) > 1 | is.infinite( a ) )
    if ( fun == "vOwenT" ) j <- j & opt & abs( a ) < 1
    pah <- h # dummy - forced to have a mpfr class if Rmpfr is used
    pah[ j ]  <- pnorm( a[ j ] * h[ j ] )
    pah[ !j ] <- 0 # h == 0 & ( a == Inf | a == -Inf )
    w <- eval( call( fun, h, a, pnorm( h ), pah, opt ) )
    z[ i ] <- w
    n[ i ] <- attr( w, "nIter" )
  }
  attr( z, "nIter" ) <- n
  return( z )
}
