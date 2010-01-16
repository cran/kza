\name{Rkzp}
\alias{Rkzp}

\title{Kolmogorov-Zurbenko Periodogram}

\description{
  Kolmogorov-Zurbenko periodogram and smoothing using DiRienzo-Zurbenko (DZ) or Neagu-Zurbenko (NZ) methods.
}

\details{
  The Kolmogorov-Zurbenko Periodogram is an estimate of the spectral density using the Kolmogorov-Zurbenko Fourier Transform (KZFT).
  This is the R scripted version that does not use FFTW (see README).
}

\usage{
Rkzp(y, m=round(length(y)/10), k=1, f=NULL)
}

\arguments{
  \item{y}{The raw data.}
  \item{m}{The width of filtering window}
  \item{k}{The number of iterations for the KZFT}
  \item{f}{The frequency that KZFT is applied at.}
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
  \code{\link{smooth_kzp}},
  \code{\link{plot.kzp}},  
}

\examples{

#example 
t<-1:5000
y<-1.1*sin(2*pi*0.3*t)+7*sin(2*pi*0.4*t)+10*rnorm(length(t))
a<-Rkzp(y,m=1000)
plot.kzp(a,c=0.03)
}
