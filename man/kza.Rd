\name{kza}
\alias{kza}
\title{Kolmogorov-Zurbenko Adaptive}
\author{Brian Close <brian.close@gmail.com> and Igor Zurbenko <IZurbenko@albany.edu>}
\description{
KZA will recover 2-dimensional or 3-dimensional image or signal buried in noise.
}
\details{
The selection of parameters of KZA depend on the nature of the data. 
This function can take a long time to run, depending on the number of dimensions
and the size of the dimensions.
}
\usage{
kza(x, m, y = NULL, k = 3, min_size = round(0.05*m), tol = 1.0e-5, impute_tails = FALSE)
}
\arguments{
  \item{x}{A vector of the time series or a matrix (2d) or an array (3d) of an image.}
  \item{m}{The window for the filter.}
  \item{y}{The filtered output from kz.}
  \item{k}{The number of iterations.}
  \item{min_size}{Minimum size of window q.}
  \item{tol}{The smallest value to accept as nonzero.}
  \item{impute_tails}{The default is to drop the tails.}
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
#######
# this is an example of detection of a break point in a time series
#######
yrs <- 20
t <- seq(0,yrs,length=yrs*365)
m <- 365

#noise
e <- rnorm(n = length(t),0,1)
trend <- seq(0,-1,length=length(t))

#signal
bkpt <- 3452
brk <- c(rep(0,bkpt),rep(0.5,length(t)-bkpt))
signal <- trend + brk

# y = seasonal + trend + break point + noise
y <- sin(2*pi*t) + signal + e

k.kz <- kz(y,m)

# kza reconstruction of the signal
k.kza <- kza(y,m,y=k.kz,min_size=10)

par(mfrow=c(2,1))
plot(y,type="l", ylim=c(-3,3))
plot(signal,type="l",ylim=c(-3,3), 
    main="Signal and KZA Reconstruction")
lines(k.kza$kza, col=4)

######################
# image detection (2d)
######################
a <- matrix(rep(0,100*100),nrow=100)
a[35:70,35:70]<-1
a <- a + matrix(rnorm(100*100,0,1),nrow=100)
y<-kz(a,m=15,k=3)
v <- kza(a,m=15,y=y,k=3,impute_tails=TRUE)

x <- seq(1,100)
y <- x
op <- par(bg = "white")

###
#noise
###
persp(x, y, a, zlab="z", zlim=c(0,2), ticktype="detailed", theta = 30, phi = 30, col = "lightblue")

###
#kza filtered
###
X11()
persp(x, y, v$kza, zlab="z", zlim=c(0,2), ticktype="detailed", theta = 30, phi = 30, col = "lightblue")

###
# another view
###
par(mfrow=c(1,2))
image(a,col=gray(seq(0,1,1/255)))
image(v$kza,col=gray(seq(0,1,1/255)))
par(mfrow=c(1,1))


}
\keyword{ts}
