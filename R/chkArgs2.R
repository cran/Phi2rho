chkArgs2 <- function( x, y, rho, opt ) {
  stopifnot( length( x ) >= 1 && length( y ) >= 1 && length( rho ) >= 1 )
  stopifnot( length( abs( rho ) <= 1 ) >= 1 )
  stopifnot( ( isa( x, "numeric" ) || isa( x, "mpfr" ) ) &&
             class( x ) == class( y ) && class( y ) == class( rho ) )
  if ( isa( x, "mpfr" ) ) {
    prec_x   <- getPrec( x )
    prec_y   <- getPrec( y )
    prec_rho <- getPrec( rho )
    stopifnot( all( c( prec_x, prec_y, prec_rho ) == prec_x[ 1 ] ) )
  }
  if ( length( x ) > 1 & length( y )   > 1 ) stopifnot( length( x ) == length( y ) )
  if ( length( x ) > 1 & length( rho ) > 1 ) stopifnot( length( x ) == length( rho ) )
  if ( length( y ) > 1 & length( rho ) > 1 ) stopifnot( length( y ) == length( rho ) )
  stopifnot( length( opt ) == 1 && ( opt == TRUE || opt == FALSE) )
  return( TRUE )
}
