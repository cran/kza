\name{kza}
\alias{kza}
\title{Kolmogorov-Zurbenko Adaptive}
\author{Brian Close <bclose@nycap.rr.com> and Igor Zurbenko <IZurbenko@albany.edu>}
\description{
KZA will recover a 2-dimensional image buried in noise.
The KZA Algorithm can recover break in a signal buried in noise.
}
\details{
The selection of parameters of KZA depend on the nature of the data.
}
\usage{
kza(x, q, k = 3, m = 0, tol = 1.0e-5)
}
\arguments{
  \item{x}{A vector of the time series or a matrix of an image.}
  \item{q}{The half length of the window size for the filter.}
  \item{k}{Number of iterations to run the filter.}
  \item{m}{Minimum size of filtering window.}
  \item{tol}{The smallest value to accept as nonzero.}
}
\references{
    I. Zurbenko, P.S. Porter, S.T. Rao, J.Y. Ku, R. Gui, R.E. Eskridge
    Detecting Discontinuities in Time Series of Upper-air Data: 
    Development and Demonstration of an Adaptive Filter Technique.
    Journal of Climate: (1996) Vol. 9, No. 12, pp. 3548 3560.
    \url{http://ams.allenpress.com/amsonline/?request=get-abstract&issn=1520-0442&volume=009&issue=12&page=3548}
    
    Kevin L. Civerolo, Elvira Brankov, S. T. Rao, Igor Zurbenko
    Assessing the impact of the acid deposition control program.
    Atmospheric Environment 35 (2001) 4135-4148
    \url{http://www.elsevier.com/locate/atmosenv}

    J.Chen, I.Zurbenko, Nonparametric Boundary detection, Communications in Statistics, 
    Theory and Methods, Vol.26, 12, 2999-3014, 1997.
}
\examples{
# image detection
a <- matrix(rep(0,100*100),nrow=100)
a[35:70,35:70]<-1
a <- a + matrix(rnorm(100*100,0,2),nrow=100)
k <- kza(a,15)

par(mfrow=c(1,2))
image(a,col=gray(seq(0,1,1/255)))
image(k,col=gray(seq(0,1,1/255)))

# another example
a <- matrix(rep(0,100*100),nrow=100)
a[35:70,35:70]<-1
a <- a + matrix(rnorm(100*100,0,2),nrow=100)
k <- kza(a,15)

par(mfrow=c(1,2))
x <- seq(1,100)
y <- x
op <- par(bg = "white")
persp(x, y, a, zlab="z", zlim=c(-4,6), ticktype="detailed", 
    theta = 30, phi = 30, col = "lightblue")
persp(x, y, k, zlab="z", zlim=c(0,2), ticktype="detailed", 
    theta = 30, phi = 30, col = "lightblue")


# this is an example of detection of a break point in a time series
# seasonal data
yrs <- 20
t <- seq(0,yrs,length=yrs*365)
q <- 365

#noise
e <- rnorm(n = length(t),0,1)
trend <- seq(0,-1,length=length(t))

#signal
bkpt <- 3452
brk <- c(rep(0,bkpt),rep(0.5,length(t)-bkpt))
signal <- trend + brk

# y = seasonal + trend + break point + noise
y <- sin(2*pi*t) + signal + e

k.kz <- kz(y,q)

# kza reconstruction of the signal
k.kza <- kza(y,q,m=10)

par(mfrow=c(1,2))
plot(y,type="l", ylim=c(-3,3))
plot(signal,type="l",ylim=c(-3,3), 
    main="Signal and KZA Reconstruction")
lines(k.kza, col=4)
}
\keyword{ts}
