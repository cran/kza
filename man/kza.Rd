\name{kza}
\alias{kza}
\title{Kolmogorov-Zurbenko Adaptive}
\description{
  Detect change points in a time series.
}
\details{
A change point is an abrupt change (change of structure due to economy, sociology, climate, legislation, etc.) 
in the smooth baseline of the original data covered by noise. The KZA Algorithm can determine the unknown time 
and value of that change. KZA can be used with missing values.
}
\usage{
kza(v, q, k = 3)
}
\arguments{
  \item{v}{A vector of the time series}
  \item{q}{The half length of the window size for the filter}
  \item{k}{Number of iterations to run the filter, default = 3}
}
\references{
    I. Zurbenko, P.S. Porter, S.T. Rao, J.Y. Ku, R. Gui, R.E. Eskridge
    Detecting Discontinuities in Time Series of Upper-air Data: Development and Demonstration of an Adaptive Filter Technique.
    Journal of Climate: (1996) Vol. 9, No. 12, pp. 3548–3560.
    \url{http://ams.allenpress.com/amsonline/?request=get-abstract&issn=1520-0442&volume=009&issue=12&page=3548}
    
    Kevin L. Civerolo, Elvira Brankov, S. T. Rao, Igor Zurbenko
    Assessing the impact of the acid deposition control program.
    Atmospheric Environment 35 (2001) 4135-4148
    \url{http://www.elsevier.com/locate/atmosenv}
}
\examples{
# period of sesonality is 365
# initial parameter of kza 365

# seasonal data
yrs <- 20
t <- seq(0,yrs,length=yrs*365)
q <- 365/2

#noise
set.seed(6); e <- rnorm(n = length(t), sd = 1.0)
trend <- seq(0,-1,length=length(t))

#signal
bkpt <- 3452
brk <- c(rep(0,bkpt),rep(0.2,length(t)-bkpt))
signal <- trend + brk

# y = seasonal + trend + break point + noise
y <- sin(2*pi*t) + signal + e

# kza reconstruction of the signal
z <- kza(y,q)
k.kza <- z[[1]]

# determine adaptive filter
f <- kzf(y,q,rep(0,365*yrs))
k.f <- f[[5]]

# plot of sample variance
s <- kzsv(k.kza,q,k.f)
k.s <- s[[1]]

par(mfrow=c(2,2))
plot(y,type="l", main="Signal + Seasonality(365) + Noise ")
plot(signal,type="l",ylim=c(-2,2), main="Signal and KZA Reconstruction (break=0.2,q=365/2)")
lines(k.kza)
plot(k.s,type="l",main="Sample Variance (break=0.2,q=365/2)")

# example with break = 0.5
brk <- c(rep(0,bkpt),rep(0.5,length(t)-bkpt))
signal <- trend + brk

# y = seasonal + trend + break point + noise
y <- sin(2*pi*t) + signal + e

# kza reconstruction of the signal
z <- kza(y,365)
k.kza <- z[[1]]
plot(signal,type="l",ylim=c(-2,2), main="Signal and KZA Reconstruction (break=0.5,q=365)")
lines(k.kza)
}
\keyword{ts}
