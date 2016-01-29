\name{kzsv}
\alias{kzsv}
\alias{summary.kzsv}
\alias{plot.kzsv}
\title{Kolmogorov-Zurbenko Adaptive filter with Sample Variance.}
\description{
  Sample variance of a Kolmogorov-Zurbenko adaptive filter.
  You want a sigma of at least 3.
}
\usage{
kzsv(object)
\method{summary}{kzsv}(object, digits=getOption("digits"), \dots)
\method{plot}{kzsv}(x, \dots)
}
\arguments{
  \item{object}{The resultant object from kza function.}
  \item{x}{The results of the kza function.}
  \item{digits}{Precision of output.}
  \item{...}{Other parameters.}
}
\examples{
x <- c(rep(0,4000),rep(0.5,2000),rep(0,4000))
noise <- rnorm(n = 10000, sd = 1.0) # normally-distributed random variates
v <- x + noise
a<-kza(v, m=1000, k=3)
sv<-kzsv(a)
}
\keyword{ts}
