\name{OwenT}
\alias{OwenT}
\title{Owen's T Function}
\description{Computes Owen's T function using the modified Euler's arctangent series, tetrachoric series or Vasicek's series.}
\usage{OwenT(h, a, opt = TRUE, fun = c("mOwenT", "tOwenT", "vOwenT"))}
\arguments{
  \item{h}{Numeric scalar or vector.}
  \item{a}{Numeric scalar or vector.}
  \item{opt}{If TRUE, an optimized calculation is performed.}
  \item{fun}{The name of the internal function being used.}
}
\details{
If \code{h} and \code{a} are both vectors, they must be of the same length. If one of \code{h} and \code{a} is a vector and the other is a scalar, the latter is replicated to the length of the former. The calculation is performed component-wise.

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
\note{Function is ready to work with the \pkg{Rmpfr} package, which enables using arbitrary precision numbers instead of double precision ones.  Assuming \pkg{Rmpfr} is loaded, it is sufficient to be called with parameters \sQuote{h} and \sQuote{a}, which have class \sQuote{mpfr} and the same precision.}
\seealso{\code{\link{Phi2xy}}}

\examples{
OwenT(2, 0.5)
OwenT(2, 0.5, fun = "mOwenT") # modified arctangent series (default)
OwenT(2, 0.5, fun = "tOwenT") # tetrachoric series
OwenT(2, 0.5, fun = "vOwenT") # Vasicek's series

rho <- 0.6
a <- rho/sqrt(1 - rho^2)
OwenT(0.3, a)
OwenT(0.3, a, fun = "tOwenT") 
OwenT(c(-1, 0.5, 4), a, fun = "vOwenT")

OwenT(2, c(-1, -0.5, 0, 0.5, 1), fun = "vOwenT")

function (h, a, opt = TRUE, fun = c("mOwenT", "tOwenT", "vOwenT")) 
{
    chkArgs1(h = h, a = a, opt = opt)
    if (length(h) < length(a)) 
        h <- rep(h, length(a))
    if (length(a) < length(h)) 
        a <- rep(a, length(h))
    z <- h
    n <- rep(0, length(z))
    z[a == 0] <- 0
    z[is.na(a)] <- NA
    i <- a != 0 & !is.na(a)
    if (any(i)) {
        h <- h[i]
        a <- a[i]
        ph <- pnorm(h)
        fun <- match.arg(fun)
        if (fun == "mOwenT") 
            j <- abs(a) > 1
        if (fun == "tOwenT") 
            j <- opt & abs(a) > 1
        if (fun == "vOwenT") 
            j <- opt & abs(a) < 1
        pah <- h
        if (any(j)) 
            pah[j] <- pnorm(a[j] * h[j])
        pah[is.nan(pah)] <- 0
        w <- eval(call(fun, h, a, ph, pah, opt))
        z[i] <- w
        n[i] <- attr(w, "nIter")
    }
    attr(z, "nIter") <- n
    return(z)
  }
}
\concept{ Owen's T function }
