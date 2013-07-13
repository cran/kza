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
kzs(y,m=NULL,k=3,t=NULL)
}
\arguments{
  \item{y}{data}
  \item{m}{smooth}
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
true<-((exp(2.5*x)+sin(25*x))-1)/3

noise <- rnorm(n)
y <- true + noise 

a<-kzs(y,m=60)

par(mfrow=c(2,1))
plot(y,type='l')
lines(true,col="red")

plot(a,type='l', ylim=c(-2,4))
lines(true,col="red")
par(mfrow=c(1,1))

################
# second example
################
t <- seq(from = -round(400*pi), to = round(400*pi), by = .25) 
ts <- 0.5*sin(sqrt((2*pi*abs(t))/200))
signal <- ifelse(t < 0, -ts, ts)
et <- rnorm(length(t), mean = 0, sd = 1)
yt <- et + signal

b<-kzs(yt,m=400)
par(mfrow=c(2,1))
plot(yt,type='l')
lines(signal,col="red")

plot(b,type='l', ylim=c(-0.5,1))
lines(signal,col="red")
par(mfrow=c(1,1))

}
