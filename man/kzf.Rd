\name{kzf}
\alias{kzf}
\title{Kolmogorov-Zurbenko adaptive filter}
\description{
  kzf determines Kolmogorov-Zurbenko adpative filter. The adaptive filter is used by kzsv for sample variance.
}
\usage{
kzf(v, q, d, k =3)
}
\arguments{
  \item{v}{A vector of the time series}
  \item{q}{The half length of the window size for the filter}
  \item{d}{adaptive filter}
  \item{k}{Number of iterations, default = 3}
}
\examples{
x <- c(rep(0,1000),rep(0.5,1000),rep(0,1000))
set.seed(1)

v <- x + rnorm(n = 3000, sd = 1.0)    # normally-distributed random variates

z <- kza(v,100)
z <- z[[1]]

d <- rep(0,3000)
f <- kzf(v,100,d)
f <- f[[5]]

s <- kzsv(z,100,f)
plot(s[[1]],type="l")
}
\keyword{ts}
