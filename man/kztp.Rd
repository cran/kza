\name{kztp}
\alias{kztp}

\title{Kolmogorov-Zurbenko Third-Order Periodogram}

\description{
  Kolmogorov-Zurbenko Third-Order Periodogram for estimating spectrums
}

\details{
  The Kolmogorov-Zurbenko Third-Order Periodogram is used to estimate spectral density of a signal. The smoothing methods are adaptive allowing the bandwidth of the spectral window to vary according to the smoothness of the underlying spectral density. For details, please see to DiRienzo and Zurbenko (1998) and Neagu and Zurbenko (2003).  
}

\usage{
kztp(x,m,k,box=c(0,0.5,0,0.5))
}

\arguments{
  \item{x}{The signal.}
  \item{m}{The window size for the kzft filter.}
  \item{k}{The number of iterations.}
  \item{box}{The window for the application of third-order periodgram.}
}

\references{
  I. G. Zurbenko, 1986: The spectral Analysis of Time Series. North-Holland, 248 pp.
  
  I. G. Zurbenko, P. S. Porter, Construction of high-resolution wavelets, Signal Processing 65: 315-327, 1998.
  
  W. Yang, I. G. Zurbenko, A semi-adaptive smoothing algorithm in bispectrum estimation, Proceedings of the American Statistical Association, Seattle, 2006.
  
  Wei Yang and Igor Zurbenko, kzft: Kolmogorov-Zurbenko Fourier Transform and Applications, R-Project 2007.
  
}

\seealso{
  \code{\link{kzft}},
  \code{\link{kzp}},
  \code{\link{smooth_kzp}},
}

\examples{
#
# example
#
t<-1:10000
y<-sin(2*pi*0.1*t)+3*sin(2*pi*0.2*t) + rnorm(length(t),0,2)

a<-kztp(y,50,1)
z<-log(Mod(a))
s<-smooth_kzp(z,0.01,20)

omega<-seq(0,1,length=51)[2:26]
#filled.contour(omega,omega,s,xlab="freq",ylab="freq",main="Smoothed 3rd Order Periodogram")
}
