\name{kzsv}
\alias{kzsv}
\title{Kolmogorov-Zurbenko Adaptive filter with Sample Variance.}
\description{
  Sample variance of a Kolmogorov-Zurbenko adaptive filter.
  You should get a sigma of at least 3 but preferably more than that.
}
\usage{
kzsv(y, kza, kz, q, k=3, m=round(0.05*q), tol=1.0e-5)
}
\arguments{
  \item{y}{The raw data.}
  \item{kza}{A vector of the resultant time series from kza function.}
  \item{kz}{A moving average result from the kz function.}
  \item{q}{The window size for the filter.}
  \item{k}{The number of iterations that were used.}
  \item{m}{Minimum size of filtering window.}
  \item{tol}{The smallest value to accept as nonzero.}
}
\examples{
x <- c(rep(0,4000),rep(0.5,2000),rep(0,4000))
noise <- rnorm(n = 10000, sd = 1.0) # normally-distributed random variates
v <- x + noise
kzsv(v, q=1000, k=3)
}
\keyword{ts}
