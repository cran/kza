\name{kz}
\alias{kz}
\title{Kolmogorov-Zurbenko filter}
\description{
  Kolmogorov-Zurbenko low-pass linear filter.
}
\details{
  KZ is an iterated moving average. The filter can be used with missing values.
}
\usage{
kz(x, q, k = 3)
}
\arguments{
  \item{x}{A vector of the time series.}
  \item{q}{Window size for the filter.}
  \item{k}{Number of iterations.}
}
\references{
  Zurbenko, I. G., 1986: The spectral Analysis of Time Series. North-Holland, 248 pp.
}
\examples{
#seperation of signals
yrs <- 20
t <- seq(0,yrs,length=yrs*365)
y <- sin(2*pi*t) + sin(3*pi*t)
k.kz <- kz(y,365/2)

par(mfrow=c(3,1))
plot(y,type="l",main="y = sin(2*pi*t)+sin(3*pi*t)")
plot(k.kz,type="l",main="KZ filter")

r <- y - 4*k.kz
plot(r,type="l",main="(y - 4*kz) ~ sin(3*pi*t)")

#another example
#remove noise and high frequency
yrs <- 20
t <- seq(0,yrs,length=yrs*365)
set.seed(6); e <- rnorm(n = length(t), sd = 1.0)
y <- sin(2*pi*t) + sin(3*pi*t) + e

k.kz <- kz(y,365/4)

par(mfrow=c(2,1))
plot(y,type="l",main="y = sin(2*pi*t)+sin(3*pi*t) + noise")
plot(k.kz,type="l",main="KZ filter")
}
\keyword{ts}
