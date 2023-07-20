chkArgs1 <-
function( h, a, opt ) {
  stopifnot( length( h ) >= 1 && length( a ) >= 1 )
  if ( length( h ) > 1 && length( a ) > 1 ) stopifnot( length( h ) == length( a ) )
  # When plotting a curve with rho on the x-axis, xlim = c( -1.1, 1.1 ), and
  # calling OwenT() with a = rho / sqrt(1 - rho^2), we do not want to abort
  # the execution but just ignore the values outside the interval [ -1, 1 ].
  # That is why we also allow NA. However, at least one component of a must
  # be correct.
  stopifnot( length( a[ !is.na( a ) ] ) >= 1 )
  a <- a[ !is.na( a ) ]
  stopifnot( ( isa( h, "numeric" ) || isa( h, "mpfr" ) ) && class( h ) == class( a ) )
  if ( isa( h, "mpfr" ) ) {
    # it is not important here, but if a contains NA and Rmpfr is loaded,
    # it changes NA to NaN and implies length( getPrec( a ) ) < length( a )
    prec_h <- getPrec( h )
    prec_a <- getPrec( a )
    stopifnot( all( c( prec_h, prec_a ) == prec_h[ 1 ] ) )
  }
  stopifnot( length( opt ) == 1 && ( opt == TRUE || opt == FALSE) )
  return( TRUE )
}
