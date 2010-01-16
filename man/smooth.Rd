\name{smooth_kzp}
\alias{smooth_kzp}

\title{Kolmogorov-Zurbenko Third-Order Periodogram Smoothing Method}

\description{
  Adaptive Smoother for Kolmogorov-Zurbenko Third-Order Periodogram for estimating spectrums
}

\details{
  The Kolmogorov-Zurbenko Third-Order Periodogram is used to estimate spectral density of a signal. The smoothing methods are adaptive allowing the bandwidth of the spectral window to vary according to the smoothness of the underlying spectral density. For details, please see to DiRienzo and Zurbenko (1998) and Neagu and Zurbenko (2003).  
}

\usage{
smooth_kzp(x, c=0.01, w=length(x), method="DZ")
}

\arguments{
  \item{x}{The third-order periodogram output from the kztp function.}
  \item{c}{The percentage of variation to be used for smoothing.}
  \item{w}{The width of the smoothing window.}
  \item{method}{The smoothing method to apply.}
}

\references{
  I. G. Zurbenko, 1986: The spectral Analysis of Time Series. North-Holland, 248 pp.
  
  I. G. Zurbenko, P. S. Porter, Construction of high-resolution wavelets, Signal Processing 65: 315-327, 1998.
  
  A. G. DiRienzo, I. G. Zurbenko, Semi-adaptive nonparametric spectral estimation, Journal of Computational and Graphical Statistics 8(1): 41-59, 1998.
  
  R. Neagu, I. G. Zurbenko, Algorithm for adaptively smoothing the log-periodgram, Journal of the Franklin Institute 340: 103-123, 2003.

  W. Yang, I. G. Zurbenko, A semi-adaptive smoothing algorithm in bispectrum estimation, Proceedings of the American Statistical Association, Seattle, 2006.
  
  Wei Yang and Igor Zurbenko, kzft: Kolmogorov-Zurbenko Fourier Transform and Applications, R-Project 2007.
  
}

\seealso{
  \code{\link{kzft}},
  \code{\link{kzp}},
  \code{\link{kztp}},
}

\examples{
#
# example
#
t<-1:5000
y<-1.1*sin(2*pi*0.3*t)+7*sin(2*pi*0.4*t)+10*rnorm(length(t))
a<-kzp(y,m=1000)
plot.kzp(a,c=0.03)
}


