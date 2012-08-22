\name{kzp}
\alias{kzp}
\alias{Rkzp}
\alias{smooth.kzp}
\alias{nonlinearity.kzp}
\alias{variation.kzp}
\alias{summary.kzp}
\alias{plot.kzp}
\title{Kolmogorov-Zurbenko Periodogram}
\description{
  Kolmogorov-Zurbenko periodogram and smoothing using DiRienzo-Zurbenko (DZ).
}
\details{
  The Kolmogorov-Zurbenko Periodogram is an estimate of the spectral density using the Kolmogorov-Zurbenko Fourier Transform (KZFT).
}
\usage{
kzp(y, m=length(y), k=1, double_frequency=FALSE)
\method{smooth}{kzp}(object, log=TRUE, smooth_level=0.05, method = "DZ")
\method{nonlinearity}{kzp}(x)
\method{variation}{kzp}(x)
\method{summary}{kzp}(object, digits=getOption("digits"), top=1, \dots)
\method{plot}{kzp}(x, \dots)
Rkzp(y, m=NULL, k=3, double_frequency=FALSE)
}
\arguments{
  \item{y}{The raw data.}
  \item{m}{The width of filtering window}
  \item{k}{The number of iterations for the KZFT}
  \item{double_frequency}{The return vector is half the width of the filtering window, setting this to true will give the second half.}
  \item{object}{Output from kzp function.}
  \item{log}{Use logarithm values for smoothing.}
  \item{smooth_level}{Percentage of smoothness to apply.}
  \item{method}{Method used for smoothing; choices are "DZ" or "NZ".}
  \item{digits}{precision of output.}
  \item{top}{list top values}
  \item{...}{Other parameters.}
  \item{x}{periodogram}
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
}
\examples{
t<-1:6000
f1<-0.03
f2<-0.04
noise<-15*rnorm(length(t))
amp=1.5
s<-amp*sin(2*pi*f1*t)+amp*sin(2*pi*f2*t)
system.time(a<-kzp(s+noise,500,k=3))
b<-smooth.kzp(a, smooth_level=0.08)
par(mfrow=c(3,1))
plot(periodogram(s+noise),type='l')
plot(a)
plot(b)
par(mfrow=c(1,1))

# signal/noise
signal<-kzft(s+noise,m=500,k=3,dim=1)
print(paste("signal-to-noise ratio = ", round(sqrt(var(2*Re(signal))/var(s+noise)),4) ))

summary(a, digits=2, top=2)
}
