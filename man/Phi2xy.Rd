\name{Phi2xy}
\alias{Phi2xy}
\title{Bivariate Normal Integral}
\description{Computes the bivariate normal integral Phi2(x, y, rho).}
\usage{Phi2xy(x, y, rho, opt = TRUE, fun = c("mOwenT", "tOwenT", "vOwenT"))}
\arguments{
  \item{x}{Numeric scalar or vector.}
  \item{y}{Numeric scalar or vector.}
  \item{rho}{Numeric scalar or vector.}
  \item{opt}{If TRUE, an optimized calculation is performed.}
  \item{fun}{The name of the internal function being used.}
}
\details{
The parameter \sQuote{rho} (or at least one of its components) must be from the interval [-1,1].

Vector parameters must be of the same length, and any scalar parameters are replicated to the same length. The calculation is performed component-wise.

The parameter \code{fun} specifies which series is used:
\describe{
   \item{\dQuote{mOwenT}:}{modified Euler's arctangent series (default).}
   \item{\dQuote{tOwenT}:}{tetrachoric series.}
   \item{\dQuote{vOwenT}:}{Vasicek's series.}
}
The \code{opt} parameter enables checking the results in the submitted article and may be dropped later.

If \code{fun = "mOwenT"} and \code{opt = TRUE}, the external arctangent function is used, otherwise all necessary values are calculated on the fly, but usually more iterations are needed.

If \code{fun = "tOwenT"} or \code{fun = "vOwenT"}, and \code{opt = TRUE}, then the parameters transformation is performed when it makes sense, which significantly reduces the number of iterations.
}
\value{The value of computed function is returned, scalar or vector.  The attribute \sQuote{nIter} of returned value means the number of iterations.}
\references{
Komelj, J. (2023): The Bivariate Normal Integral via Owen's T Function as a Modified Euler's Arctangent Series, \emph{American Journal of Computational Mathematics}, \bold{13}, 4, 476--504, \doi{10.4236/ajcm.2023.134026} (or reprint \url{https://arxiv.org/pdf/2312.00011.pdf} with better typography).

Owen, D. B. (1956): Tables for Computing Bivariate Normal Probabilities, \emph{The Annals of Mathematical Statistics}, \bold{27}, 4, 1075--1090, \doi{10.1214/aoms/1177728074}.

Owen, D. B. (1980): A table of normal integrals, \emph{Communications in Statistics -- Simulation and Computation}, \bold{9}, 4, 389--419, \doi{10.1080/03610918008812164}.
}

\author{Janez Komelj}
\note{Function is ready to work with the \pkg{Rmpfr} package, which enables using arbitrary precision numbers instead of double precision ones.  Assuming \pkg{Rmpfr} is loaded, it is sufficient to be called with parameters \sQuote{x}, \sQuote{y} and \sQuote{rho}, which have class \sQuote{mpfr} and the same precision.}
\seealso{\code{\link{OwenT}}}

\examples{
Phi2xy(2, 1.3, 0.5) 
Phi2xy(-2, 0.5, -0.3, fun = "tOwenT")
Phi2xy(c(1, 2, -1.5), c(-1, 1, 2.3), 0.5, fun = "vOwenT")
Phi2xy(1, 2, c(-1, -0.5, 0, 0.5, 1))
Phi2xy(c(1, 2), c(-1,3), c(-0.5, 0.8))

function (x, y, rho, opt = TRUE, fun = c("mOwenT", "tOwenT", 
    "vOwenT")) 
{
    chkArgs2(x = x, y = y, rho = rho, opt = opt)
    fun <- match.arg(fun)
    sgn <- function(x) {
        y <- sign(x)
        y[y == 0] <- 1
        return(y)
    }
    frx <- function(x, y, rho) {
        rx <- (y - rho * x)/(x * sqrt(1 - rho^2))
        rx[is.nan(rx)] <- 0
        return(-abs(rx) * sgn(y - rho * x) * sgn(x))
    }
    fprx <- function(r, x) {
        u <- r * x
        j <- !is.nan(u)
        if (fun == "mOwenT") 
            j <- j & abs(r) > 1
        if (fun == "tOwenT") 
            j <- j & opt & abs(r) > 1
        if (fun == "vOwenT") 
            j <- j & opt & abs(r) < 1
        u[j] <- pnorm(u[j])
        u[is.nan(u)] <- 0
        return(u)
    }
    fz <- function(x, y, rho, rx, ry, px, py) {
        n <- rep(0, length(x))
        i <- rho != 0 & abs(rho) < 1 & (x != 0 | y != 0)
        z <- x
        if (any(!i)) {
            z[rho == 0] <- px[rho == 0] * py[rho == 0]
            z[rho == +1] <- pmin(px[rho == +1], py[rho == +1])
            z[rho == -1] <- pmax(px[rho == -1] + py[rho == -1] - 
                1, 0)
            j <- rho != 0 & abs(rho) < 1
            if (any(j)) 
                z[j] <- 1/4 + asin(rho[j])/(2 * pi)
        }
        if (any(i)) {
            prx <- fprx(rx[i], x[i])
            pry <- fprx(ry[i], y[i])
            zz <- eval(call(fun, c(x[i], y[i]), c(rx[i], ry[i]), 
                c(px[i], py[i]), c(prx, pry), opt, TRUE))
            n[i] <- attr(zz, "nIter")
            z[i] <- (px[i] + py[i])/2 + zz
            j <- i & (x * y < 0 | x * y == 0 & x + y < 0)
            z[j] <- z[j] - 1/2
        }
        attr(z, "nIter") <- n
        return(z)
    }
    dim <- max(length(x), length(y), length(rho))
    if (isa(x, "mpfr")) {
        pi <- Rmpfr::Const("pi", Rmpfr::getPrec(x))
        z <- mpfrArray(NA, dim = dim, precBits = Rmpfr::getPrec(x))
    }
    else z <- array(NA, dim = dim)
    n <- array(0, dim = dim)
    px <- pnorm(x)
    py <- pnorm(y)
    if (length(x) < dim) 
        x <- rep(x, dim)
    if (length(y) < dim) 
        y <- rep(y, dim)
    if (length(rho) < dim) 
        rho <- rep(rho, dim)
    if (length(px) < dim) 
        px <- rep(px, dim)
    if (length(py) < dim) 
        py <- rep(py, dim)
    k <- !is.na(rho) & abs(rho) <= 1
    x <- x[k]
    y <- y[k]
    rho <- rho[k]
    px <- px[k]
    py <- py[k]
    q <- (x^2 - 2 * rho * x * y + y^2)/(2 * (1 - rho^2))
    phi <- exp(-q)/(2 * pi * sqrt(1 - rho^2))
    phi[is.nan(phi)] <- 0
    i <- phi > 1
    j <- rho[i] < 0
    n1 <- length(rho[i])
    n2 <- length(rho) - n1
    dim <- 2 * n1 + n2
    if (isa(x, "mpfr")) 
        xx <- mpfrArray(0, dim = dim, precBits = Rmpfr::getPrec(x))
    else xx <- rep(0, dim)
    yy <- xx
    rr <- xx
    pxx <- xx
    pyy <- xx
    if (n1 > 0) {
        r <- rho[i]
        u <- x[i]
        v <- y[i] * sgn(r)
        w <- (u - v)/sqrt(2 * (1 - abs(r)))
        r <- -sqrt((1 - abs(r))/2)
        i1 <- 1:(2 * n1)
        xx[i1] <- c(w, -w)
        yy[i1] <- c(v, u)
        rr[i1] <- c(r, r)
        pw <- pnorm(w)
        pv <- (1 - sgn(rho[i]))/2 + sgn(rho[i]) * py[i]
        pxx[i1] <- c(pw, 1 - pw)
        pyy[i1] <- c(pv, px[i])
    }
    if (n2 > 0) {
        i2 <- (2 * n1 + 1):(2 * n1 + n2)
        xx[i2] <- x[!i]
        yy[i2] <- y[!i]
        rr[i2] <- rho[!i]
        pxx[i2] <- px[!i]
        pyy[i2] <- py[!i]
    }
    rx <- frx(xx, yy, rr)
    ry <- frx(yy, xx, rr)
    zz <- fz(xx, yy, rr, rx, ry, pxx, pyy)
    nn <- attr(zz, "nIter")
    if (n1 > 0) {
        s <- zz[1:n1] + zz[(n1 + 1):(2 * n1)]
        s[j] <- px[i][j] - s[j]
        z[k][i] <- s
        n[k][i] <- nn[1:n1] + nn[(n1 + 1):(2 * n1)]
    }
    if (n2 > 0) {
        z[k][!i] <- zz[i2]
        n[k][!i] <- nn[i2]
    }
    attr(z, "nIter") <- n
    return(z)
  }
}
\concept{ bivariate normal integral }
