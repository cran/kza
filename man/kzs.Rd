\name{kzs}
\alias{kzs}
\title{Kolmogorov-Zurbenko Spline}
\description{
Kolmogorov-Zurbenko Spline
}
\details{
Kolmogorov-Zurbenko Spline is essentially the Kolmogorov-Zurbenko Fourier Transform at the zero frequency.
}
\usage{
kzs(y, m, k = 1, t = NULL)

}
\arguments{
  \item{y}{The raw time series.}
  \item{m}{The window size for transform}
  \item{k}{The number of iterations for applying the KZFT}
  \item{t}{An indexing set}
}

\references{
  I. G. Zurbenko, The spectral Analysis of Time Series. North-Holland, 1986.
  
  I. G. Zurbenko, P. S. Porter, Construction of high-resolution wavelets, Signal Processing 65: 315-327, 1998.
  
  R. H. Shumway, D. S. Stoffer, Time Series Analysis and Its Applications: With R Examples, Springer, 2006.
  
  Derek Cyr and Igor Zurbenko, kzs: Kolmogorov-Zurbenko Spatial Smoothing and Applications, R-Project 2008.
  
  Derek Cyr and Igor Zurbenko, A Spatial Spline Algorithm and an Application to Climate Waves Over the United States, 2008 Joint Statistical Meetings.
}

\seealso{
  \code{\link{kzft}},
}

\examples{
n <- 1000
x <- (1:n)/n
true<-((exp(1.5*x)+sin(7*x))-1)/3

noise <- rnorm(n)
y <- true + noise 

a<-kzs(y,m=100)

plot(true,type='l')
lines(a,col="blue")
}
