\name{kzsv}
\alias{kzsv}
\title{Kolmogorov-Zurbenko Adaptive filter with Sample Variance.}
\description{
  Sample variance of a Kolmogorov-Zurbenko adaptive filter.
}
\usage{
kzsv(v, q, f)
}
\arguments{
  \item{v}{A vector of the resultant time series from kza function.}
  \item{q}{The half length of the window size for the filter.}
  \item{f}{adaptive filter derived from kzf function.}
}
\examples{
x <- c(rep(0,1000),rep(0.5,1000),rep(0,1000))
set.seed(1)

v <- x + rnorm(n = 3000, sd = 1.0)    # normally-distributed random variates

z <- kza(v,100,3)
z <- z[[1]]

# determine adaptive filter
d <- rep(0,3000)
f <- kzf(v,100,d)
f <- f[[4]]

# plot of sample variance
s <- kzsv(z,100,f)
plot(s[[1]],type="l")
}
\keyword{ts}
