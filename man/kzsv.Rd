\name{kzsv}
\alias{kzsv}
\title{Kolmogorov-Zurbenko Adaptive filter with Sample Variance.}
\description{
  Sample variance of a Kolmogorov-Zurbenko adaptive filter.
}
\usage{
kzsv(y=NULL, kza=NULL, kz=NULL, q, k=3, m=0, tol=1.0e-5)
}
\arguments{
  \item{y}{The raw data to be analyzed.}
  \item{kza}{A vector of the resultant time series from kza function.}
  \item{kz}{A moving average result from the kz function.}
  \item{q}{The half length of the window size for the filter.}
  \item{k}{The number of iterations to run kz and kza functions.}
  \item{m}{Minimum size of filtering window.}
  \item{tol}{The smallest value to accept as nonzero.}
}
\examples{
x <- c(rep(0,1000),rep(0.5,1000),rep(0,1000))
v <- x + rnorm(n = 3000, sd = 1.0)    # normally-distributed random variates
z <- kza(v,q=200,m=10)
par(mfrow=c(3,1))
plot(v,type="l")
plot(z,type="l")

s <- kzsv(y=z,kza=z,q=200)
plot(s/mean(s),type="l")
}
\keyword{ts}
