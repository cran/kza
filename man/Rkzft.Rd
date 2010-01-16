\name{Rkzft}
\alias{Rkzft}
\title{Kolmogorov-Zurbenko Fourier Transform}
\description{
Kolmogorov-Zurbenko Fourier Transform is an iterated Fourier transform.
}
\details{
Kolmogorov-Zurbenko Fourier Transform (KZFT) is the Fourier transform applied over every segment of length m iterated k times. If a frequency is not selected then the frequency with the largest amplitude is used.
The approach taken here is to iterate fft to create a matrix of frequencies over time, another approach is to use matrices for the KZFT transform (see the KZFT package by Wei Yang in R).

This is the R scripted version of kzft, it does not use FFTW (see README).

}
\usage{
Rkzft(x, m=100, k=1, f=NULL, dim=1)
}
\arguments{
  \item{x}{The raw time series}
  \item{m}{The window size for transform}
  \item{k}{The number of iterations for applying the KZFT}
  \item{f}{The frequency that KZFT is applied at.}
  \item{dim}{A value of 1 will return a vector of the given frequency and a value of 2 will return a matrix (spectra).}
}

\references{
  I. G. Zurbenko, The spectral Analysis of Time Series. North-Holland, 1986.
  
  I. G. Zurbenko, P. S. Porter, Construction of high-resolution wavelets, Signal Processing 65: 315-327, 1998.
  
  R. Neagu, I. G. Zurbenko, Tracking and separating non-stationary multi-component chirp signals, J. Franklin Inst., 499-520, 2002.
  
  R. H. Shumway, D. S. Stoffer, Time Series Analysis and Its Applications: With R Examples, Springer, 2006.
  
  Wei Yang and Igor Zurbenko, kzft: Kolmogorov-Zurbenko Fourier Transform and Applications, R-Project 2007.
  
}

\seealso{
  \code{\link{kzp}},
  \code{\link{kztp}},
}

\examples{
# example taken from Wei Yang's KZFT package
# coefficients of kzft(201,5)


# function to calculate polnomial coefficients for kzft
#require(polynom)
#coeff <- function(m, k)
#{
#    poly<-polynomial(rep(1,m))
#    polyk<-poly^k
#    coef<-as.vector(polyk)
#    coef<-coef/m^k
#    M=(m-1)*k+1
#    return(coef[1:M])
#}

#a<-coeff(201,5);
#t<-seq(1:1001)-501;
#z<-cos(2*pi*0.025*t);
#plot(z*a,type="l",xlab="Time", ylab="Coefficient", main="Coefficients of the kzft");
#lines(a);
#lines(-1*a);

# example taken from Wei Yang's KZFT package
# transfer function of the kzft(201,5) at frequency 0.025
lamda<-seq(-0.1,0.1,by=0.001)
tf1<-transfer_function(201,1,lamda,0.025)
tf2<-transfer_function(201,5,lamda,0.025)
matplot(lamda,cbind(log(tf1),log(tf2)),type="l",ylim=c(-15,0),
	ylab="Natural log transformation of the coefficients", 
	xlab="Frequency (cycles/time unit)",
    main="Transfer function of Rkzft(201,5) at frequency 0.025")

# example 
# signal with a frequency of 0.01
t<-1:1000
f<-10/1000
x<-cos(2*pi*f*t) + rnorm(length(t),0,3)

system.time(z1<-Rkzft(x,m=200,k=1, dim=1))
system.time(z2<-Rkzft(x,m=200,k=2, dim=1))
system.time(z3<-Rkzft(x,m=200,k=3, dim=1))

par(mfrow=c(2,2))
plot(x,type="l",main="Original time series",xlab="t", ylab="y")
plot(2*Re(z1),type="l",main="kzft(200,1)",xlab="t", ylab="y")
plot(2*Re(z2),type="l",main="kzft(200,2)",xlab="t", ylab="y")
plot(2*Re(z3),type="l",main="kzft(200,3)",xlab="t", ylab="y")
par(mfrow=c(1,1))
}
