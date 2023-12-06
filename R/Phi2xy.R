Phi2xy <- function( x, y, rho, opt = TRUE,
                    fun = c( "mOwenT", "tOwenT", "vOwenT" ) ) {
  # The function computes the bivariate normal integral Phi_2( x, y, rho ).
  # Any vector parameters must be of the same length. The scalar ones are
  # replicated to this length and a component-wise calculation is performed.
  # If opt == TRUE, the transformed (accelerated) tetrachoric and Vasicek''s
  # series are used if fun == "tOwenT" and fun == "vOwenT" respectively, and
  # an external arctangent function is used if fun == "mOwenT".
  chkArgs2( x = x, y = y, rho = rho, opt = opt )
  fun <- match.arg( fun )
  sgn <- function( x ) { y <- sign( x ); y[ y == 0 ] <- 1; return( y ) }
  frx <- function( x, y, rho ) { # rx calculation (rx is Owen''s ax)
    # all( abs( rho ) <= 1 ) is TRUE here
    # Attention: if abs( rho ) = 1 and x < 0, the sign of x in denominator is
    # lost and we get Inf with the wrong sign. If sign( x ) is used, we can
    # get NaN (Inf * 0) as a result instead of Inf because sign( 0 ) = 0.   
    rx <- ( y - rho * x) / ( x * sqrt( 1 - rho^2 ) )
    # Meyer <- TRUE
    # if ( Meyer ) {
    #   r  <- ( x - y * sgn( rho ) ) / ( x * sqrt( 1 - rho^2 ) )
    #   rx <- r - sqrt( ( 1 - abs( rho ) ) / ( 1 + abs( rho ) ) )     
    # }
    rx[ is.nan( rx ) ] <- 0 # 0 / 0 if x and y are 0
    # The minus sign is changed to plus in the Owen''s original formula
    # Phi_2( x, y, rho ) = ( Phi( x ) + Phi( y ) ) / 2 - T( x, rx ) - T( y, rx )
    # and since T is an odd function of the second parameter, the sign of rx is
    # also changed below.
    return( -abs( rx ) * sgn( y - rho * x ) * sgn( x ) )
  }
  fprx <- function( r, x ) { # pnorm( r * x ) calculation
    u <- r * x # can be Inf * 0 or -Inf * 0
    j <- !is.nan( u )
    if ( fun == "mOwenT" ) j <- j & abs( r ) > 1
    if ( fun == "tOwenT" ) j <- j & ( opt & abs( r ) > 1 | is.infinite( r ) )
    if ( fun == "vOwenT" ) j <- j & opt & abs( r ) < 1
    u[ j ]  <- pnorm( u[ j ] )
    u[ !j ] <- 0  # ( r == Inf | r == -Inf ) & x == 0
    return( u )
  }
  fz <- function( x, y, rho, rx, ry, px, py ) {
    # all parameters are of the same length and all( abs( rho ) <= 1 ) is TRUE
    n <- rep( 0, length( x ) )
    i <- rho != 0 & abs( rho ) < 1 & ( x != 0 | y != 0 )
    z <- x
    if ( any( !i ) ) { # special cases
      z[ rho ==  0 ] <- px[ rho == 0 ] * py[ rho == 0 ]
      z[ rho == +1 ] <- pmin( px[ rho == +1 ], py[ rho == +1 ] )
      z[ rho == -1 ] <- pmax( px[ rho == -1 ] + py[ rho == -1 ] - 1, 0 )
      j <- rho != 0 & abs( rho ) < 1 # x == 0 & y == 0
      if ( any( j ) ) z[ j ] <- 1 / 4 + asin( rho[ j ] ) / ( 2 * pi )
    }
    if ( any( i ) ) { # ordinary cases
      prx <- fprx( rx[ i ], x[ i ] )
      pry <- fprx( ry[ i ], y[ i ] )
      # T( x, rx ) + T( y, ry ) computation
      zz <- eval( call( fun, c( x[ i ], y[ i ] ), c( rx[ i ], ry[ i ] ), 
                        c( px[ i ], py[ i ] ), c( prx, pry ), opt, TRUE ) )
      n[ i ] <- attr( zz, "nIter" )
      z[ i ] <- ( px[ i ] + py[ i ] ) / 2 + zz # - zz in the Owen''s formula
      j <- i & ( x * y < 0 | x * y == 0 & x + y < 0 )
      z[ j ] <- z[ j ] - 1 / 2
    }
    attr( z, "nIter" ) <- n
    return( z )
  }
  dim <- max( length( x ), length( y ), length( rho ) )
  if ( isa( x, "mpfr" ) ) {
    pi <- Rmpfr::Const( "pi", Rmpfr::getPrec( x ) )
    z  <- mpfrArray( NA, dim = dim, precBits = Rmpfr::getPrec( x ) )
  } else z <- array( NA, dim = dim )
  n  <- array( 0, dim = dim ) # frame for the number of iterations
  px <- pnorm( x )
  py <- pnorm( y )
  if ( length( x )   < dim ) x   <- rep( x,   dim )
  if ( length( y )   < dim ) y   <- rep( y,   dim )
  if ( length( rho ) < dim ) rho <- rep( rho, dim )
  if ( length( px )  < dim ) px  <- rep( px,  dim )
  if ( length( py )  < dim ) py  <- rep( py,  dim )
  k <- !is.na( rho ) & abs( rho ) <= 1
  x <- x[ k ]
  y <- y[ k ]
  rho <- rho[ k ]
  px  <- px[ k ]
  py  <- py[ k ]
  q   <- ( x^2 - 2 * rho * x * y + y^2 ) / ( 2 * ( 1 - rho^2 ) )
  phi <- exp( -q ) / ( 2 * pi * sqrt( 1 - rho^2 ) )    
  phi[ is.nan( phi ) ] <- 0 # irrelevant if abs(rho) == 1
  i   <- phi > 1
  j   <- rho[ i ] < 0
  n1  <- length( rho[ i ] )
  n2  <- length( rho ) - n1
  dim <- 2 * n1 + n2
  if ( isa( x, "mpfr" ) )
    xx <- mpfrArray( 0, dim = dim, precBits = Rmpfr::getPrec( x ) )
  else
    xx <- rep( 0, dim )
  yy  <- xx
  rr  <- xx
  pxx <- xx
  pyy <- xx
  if ( n1 > 0 ) { # potentially critical cases are split into two non-critical
    r <- rho[ i ]
    u <- x[ i ]
    v <- y[ i ] * sgn( r )
    w <- ( u - v ) / sqrt( 2 * ( 1 - abs( r ) ) )
    r <- -sqrt( ( 1 - abs( r ) ) / 2 )
    i1 <- 1 : ( 2 * n1 )
    xx[ i1 ] <- c( w, -w )
    yy[ i1 ] <- c( v, u )
    rr[ i1 ] <- c( r, r )
    pw <- pnorm( w )
    pv <- ( 1 - sgn( rho[ i ] ) ) / 2 + sgn( rho[ i ] ) * py[ i ]
    pxx[ i1 ] <- c( pw, 1 - pw )
    pyy[ i1 ] <- c( pv, px[ i ] )
  }
  if ( n2 > 0 ) {
    i2 <- ( 2 * n1 + 1 ) : ( 2 * n1 + n2 )
    xx[ i2 ]  <- x[ !i ]
    yy[ i2 ]  <- y[ !i ]
    rr[ i2 ]  <- rho[ !i ]
    pxx[ i2 ] <- px[ !i ]
    pyy[ i2 ] <- py[ !i ]
  }
  rx <- frx( xx, yy, rr )
  ry <- frx( yy, xx, rr )
  zz <- fz( xx, yy, rr, rx, ry, pxx, pyy )
  nn <- attr( zz, "nIter" )
  if ( n1 > 0 ) {
    s <- zz[ 1 : n1 ] + zz[ ( n1 + 1 ) : ( 2 * n1 ) ]
    s[ j ] <- px[ i ][ j ] - s[ j ]
    z[ k ][ i ] <- s
    n[ k ][ i ] <- nn[ 1 : n1 ] + nn[ ( n1 + 1 ) : ( 2 * n1 ) ]
  }
  if ( n2 > 0 ) {
    z[ k ][ !i ] <- zz[ i2 ]
    n[ k ][ !i ] <- nn[ i2 ]
  }
  attr( z, "nIter" ) <- n
  return( z )
}
