\name{kzp}
\alias{kzp}

\title{Kolmogorov-Zurbenko Periodogram}

\description{
  Kolmogorov-Zurbenko periodogram and smoothing using DiRienzo-Zurbenko (DZ) or Neagu-Zurbenko (NZ) methods.
}

\details{
  The Kolmogorov-Zurbenko Periodogram is an estimate of the spectral density using the Kolmogorov-Zurbenko Fourier Transform (KZFT).
}

\usage{
kzp(z, m, k)
}

\arguments{
  \item{z}{The matrix output by KZFT}
  \item{m}{The width of window}
  \item{k}{The number of iterations for the KZFT}
}

\references{
  I. G. Zurbenko, 1986: The spectral Analysis of Time Series. North-Holland, 248 pp.
  
  I. G. Zurbenko, P. S. Porter, Construction of high-resolution wavelets, Signal Processing 65: 315-327, 1998.
  
  A. G. DiRienzo, I. G. Zurbenko, Semi-adaptive nonparametric spectral estimation, Journal of Computational and Graphical Statistics 8(1): 41-59, 1998.
  
  R. Neagu, I. G. Zurbenko, Algorithm for adaptively smoothing the log-periodgram, Journal of the Franklin Institute 340: 103-123, 2003.
  
  Wei Yang and Igor Zurbenko, kzft: Kolmogorov-Zurbenko Fourier Transform and Applications, R-Project 2007.
  
}

\seealso{
  \code{\link{kzft}},
  \code{\link{kztp}},
  \code{\link{smooth}},
}

\examples{

#example 
t<-1:5000
y<-1.1*sin(2*pi*0.3*t)+7*sin(2*pi*0.4*t)+10*rnorm(length(t))
m=1000
z <- kzft(y, m, 1, dim=2)
a<-log(kzp(z$Complex,m,1))
spg.dz<-smooth(a,0.03, method="DZ")
spg.nz<-smooth(a,0.03, method="NZ")

omega<-seq(0,1,length=1000)[1:500]
par(mfrow=c(2,2))
plot(omega,a,main="Raw periodogram",type="l", xlab="Frequency", ylab="")
plot(omega,spg.dz,main="Smoothed Periodogram DZ method",type="l", xlab="Frequency", ylab="")
plot(omega,spg.nz,main="Smoothed Periodogram NZ method",type="l", xlab="Frequency", ylab="")

}
