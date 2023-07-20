OwenT <-
function( h, a, opt = TRUE, fun = c( "mOwenT", "tOwenT", "vOwenT" ) ) {
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
    ph  <- pnorm( h )  
    fun <- match.arg( fun ) 
    if ( fun == "mOwenT" ) j <- abs( a ) > 1       # opt means atanExt in mOwenT
    if ( fun == "tOwenT" ) j <- opt & abs( a ) > 1 # opt means transf in tOwenT
    if ( fun == "vOwenT" ) j <- opt & abs( a ) < 1 # opt means transf in vOwenT
    pah <- h # dummy - forced to have a mpfr class if Rmpfr is used
    if ( any( j ) ) pah[ j ] <- pnorm( a[ j ] * h[ j ] )
    pah[ is.nan( pah ) ] <- 0 # ( a == Inf | a == -Inf ) & h == 0
    w <- eval( call( fun, h, a, ph, pah, opt ) )
    z[ i ] <- w
    n[ i ] <- attr( w, "nIter" )
  }
  attr( z, "nIter" ) <- n
  return( z )
}
