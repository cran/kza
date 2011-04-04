\name{periodogram}
\alias{periodogram}
\title{Periodogram}
\description{
  Raw periodogram.
}
\details{
  Periodogram is an estimate of the spectral density using FFT.
}
\usage{
periodogram(y)
}
\arguments{
  \item{y}{The raw data.}
}
\seealso{
  \code{\link{kzp}},
}
\examples{
t<-1:1000
f1<-0.3
f2<-0.4
noise<-15*rnorm(length(t))
s<-3*sin(2*pi*f1*t)+3*sin(2*pi*f2*t)
plot(periodogram(s+noise),type='l')
}
