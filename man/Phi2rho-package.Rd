\name{Phi2rho-package}
\alias{Phi2rho-package}
\alias{Phi2rho}
\docType{package}
\title{Owen's T Function and Bivariate Normal Integral}
\description{Computes the Owen's T function or the bivariate normal integral.}
\details{
The DESCRIPTION file:

Package: Phi2rho\cr
Type: Package\cr
Title: Owen's T Function and Bivariate Normal Integral\cr
Version: 1.0\cr
Date: 2023-07-19\cr
Authors@R: person("Janez","Komelj", role = c("aut","cre"),\cr
                  email = "jkomelj@siol.net")\cr
Depends: R (>= 3.5.0), stats\cr
Imports: Rmpfr\cr
Description: Computes the Owen's T function or the bivariate\cr
             normal integral using one of the following:\cr
             modified Euler's arctangent series, tetrachoric\cr
             series, or Vasicek's series.\cr
License: GPL-2 | GPL-3\cr
}
\author{Janez Komelj

Maintainer: Janez Komelj <jkomelj@siol.net>}
\references{
Komelj J.: The bivariate normal integral via Owen's T function as a modified Euler's arctangent series, submitted to ...

Owen, D. B. (1956): Tables for Computing Bivariate Normal Probabilities, \emph{The Annals of Mathematical Statistics}, \bold{27}, 4, 1075--1090, \doi{10.1214/aoms/1177728074}.

Owen, D. B. (1980): A table of normal integrals, \emph{Communications in Statistics - Simulation and Computation}, \bold{9}, 4, 389--419, \doi{10.1080/03610918008812164}.

Vasicek, O. A. (1998): A series expansion for the bivariate normal integral, \emph{The Journal of Computational Finance}, \bold{1}, 4, 5--10, \doi{10.21314/JCF.1998.015}.
}

\concept{ Owen's T function }
\concept{ bivariate normal integral }
\concept{ tetrachoric series }
\concept{ Vasicek's series }
\concept{ package }
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

Phi2xy(2, 1.3, 0.5) 
Phi2xy(-2, 0.5, -0.3, fun = "tOwenT")
Phi2xy(c(1, 2, -1.5), c(-1, 1, 2.3), 0.5, fun = "vOwenT")
Phi2xy(1, 2, c(-1, -0.5, 0, 0.5, 1))
Phi2xy(c(1, 2), c(-1,3), c(-0.5, 0.8))
}
