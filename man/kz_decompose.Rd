\name{kz_decompose}
\alias{kz_decompose}
\alias{fitted.kzd}
\alias{residuals.kzd}
\alias{trend.kzd}
\alias{frequency.kzd}
\alias{predict.kzd}
\alias{plot.kzd}
\alias{summary.kzd}
\title{Kolmogorov-Zurbenko decompose}
\description{
  Kolmogorov-Zurbenko decomposition of univariate time series into long term and seasonal components.
}
\details{
  KZ is an iterated moving average filter. 
}
\usage{
kz_decompose(y, k=3, period=NULL)
\method{fitted}{kzd}(object, \dots)
\method{residuals}{kzd}(object, \dots)
\method{trend}{kzd}(object, \dots)
\method{frequency}{kzd}(x, \dots)
\method{predict}{kzd}(object, n.ahead = 1, prediction.interval = FALSE, level = 0.95, \dots)
\method{plot}{kzd}(x, \dots)
\method{summary}{kzd}(object, digits=getOption("digits"), \dots)
}
\arguments{
  \item{y}{univariate time series.}
  \item{k}{Number of iterations.}
  \item{period}{Vector with the window sizes for the long term and seasonal period.}
  \item{object}{kzd object returned from kz_decompose.}
  \item{x}{kzd object returned from kz_decompose.}
  \item{...}{Other parameters.}
  \item{n.ahead}{The number of data points ahead to predict.}
  \item{prediction.interval}{prediction interval level}
  \item{level}{level}
  \item{digits}{precision}
}
\examples{
yrs <- 20
t <- seq(0,yrs,length=yrs*365)
amp<-10
y <- amp*sin(2*pi*t*12+10)
y<-ts(y,start=c(1980,10),frequency=365)
signal<-y

tr<-seq(1, amp/2, length.out=7300)
y.orig<-signal+exp(tr)+300
y<-y.orig
set.seed(6); e <- rnorm(n = length(t), sd = 50)
SNR<-mean(abs(y)^2)/(mean((abs(e)^2)))
y<-log(y+e)

a<-kz_decompose(y, period=c(365,30))
plot(a)
}
\keyword{kz_decompose}
