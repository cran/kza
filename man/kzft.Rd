\name{kzft}
\alias{kzft}
\alias{coeff}
\alias{max_freq}
\alias{transfer_function}
\title{Kolmogorov-Zurbenko Fourier Transform}
\description{
Kolmogorov-Zurbenko Fourier Transform is an iterated Fourier transform.
}
\details{
Kolmogorov-Zurbenko Fourier Transform (KZFT) is the Fourier transform applied 
over every segment of length m iterated k times. The argument alg="F" will use 
Fast Fourier Transforms written in C (fftw library). The alg="C" is a slow 
Fourier Transform but has the advantage of being able to handle missing values. 
It currently works in one dimension. The alg="R" is an R version of KZFT for 
experimental purposes.The coeff function generates the coefficients for the KZFT 
function.

You will introduce a phase shift and decrease the fidelity of the signal if the product of f*m is not an integer.
}
\usage{
kzft(x, m = NULL, k = 1, f = NULL, dim = NULL, index = NULL, alg = c("F","C","R"))
coeff(m, k)
max_freq(x, m, start = 1, t = NULL)
transfer_function(m, k, lamda = seq(-0.5,0.5,by=0.01), omega = 0 )
}
\arguments{
  \item{x}{The raw time series}
  \item{m}{The window size for transform}
  \item{k}{The number of iterations for applying the KZFT}
  \item{f}{The frequency that KZFT is applied at.}
  \item{dim}{A value of 1 will return a vector of the given frequency and a value of 2 will return a matrix (spectra).}
  \item{index}{An indexing set}
	\item{alg}{an option to choose different algorithms
    \itemize{
       \item \code{"C"} - a version written in C that uses a slow Fourier
       Transform, but has the advantage of handling missing values.
       \item \code{"F"} - this is written in C and uses Fast Fourier Transforms.
       \item \code{"R"} - much slower code written in R. Useful for 
         debugging and as documentation.
     }
  }
  \item{t}{An indexing set}
  \item{lamda}{The frequencies used for the calculating the transfer function.}
  \item{omega}{The frequency that KZFT is applied at.}
  \item{start}{The starting values for the index.}
}
\references{
  I. G. Zurbenko, The spectral Analysis of Time Series. North-Holland, 1986.
  
  I. G. Zurbenko, P. S. Porter, Construction of high-resolution wavelets, Signal Processing 65: 315-327, 1998.
  
  R. Neagu, I. G. Zurbenko, Tracking and separating non-stationary multi-component chirp signals, J. Franklin Inst., 499-520, 2002.
  
  R. H. Shumway, D. S. Stoffer, Time Series Analysis and Its Applications: With R Examples, Springer, 2006.
  
  Wei Yang and Igor Zurbenko, kzft: Kolmogorov-Zurbenko Fourier Transform and Applications, R-Project 2007.
  
  Igor G. Zurbenko, Amy L. Potrzeba, Tidal Waves in Atmosphere and Their Effects, Acta Geophysica Volume 58, Number 2, 356-373
  
}
\seealso{
  \code{\link{kzp}},
  \code{\link{kztp}},
}
\examples{
# example taken from Wei Yang's KZFT package
# coefficients of kzft(201,5)

# function to calculate polnomial coefficients for kzft
\dontrun{
a<-coeff(201,5);
t<-seq(1:1001)-501;
z<-cos(2*pi*0.025*t);
plot(z*a,type="l",xlab="Time", ylab="Coefficient", main="Coefficients of the kzft");
lines(a);
lines(-1*a);
}

# example taken from Wei Yang's KZFT package
# transfer function of the kzft(201,5) at frequency 0.025
lamda<-seq(-0.1,0.1,by=0.001)
tf1<-transfer_function(201,1,lamda,0.025)
tf2<-transfer_function(201,5,lamda,0.025)
matplot(lamda,cbind(log(tf1),log(tf2)),type="l",ylim=c(-15,0),
	ylab="Natural log transformation of the coefficients", 
	xlab="Frequency (cycles/time unit)",
    main="Transfer function of kzft(201,5) at frequency 0.025")

# example with missing values
set.seed(2)
period=101
f<-1/period
t<-1:2000
s<-1*sin(2*pi*f*t)
x<-s
noise<-3*rnorm(length(t))
x<-s+noise
m=101

rand_idx <- sample(t,100,replace=FALSE)
x[rand_idx]<-NA
t[rand_idx]<-NA
x<-as.vector(na.omit(x))
t<-as.vector(na.omit(t))

system.time(z1<-kzft(x, m=m, k=1, f=f, dim=1, index=t, alg="C"))
system.time(z2<-kzft(x, m=m, k=2, f=f, dim=1, index=t, alg="C"))
system.time(z3<-kzft(x, m=m, k=3, f=f, dim=1, index=t, alg="C"))

par(mfrow=c(2,2))
plot(x,type="l",main="Original time series",xlab="t", ylab="y")
lines(s,col="blue")
plot(2*Re(z1),type="l",main="kzft(101,1)",xlab="t", ylab="y", ylim=c(-6,6))
lines(s,col="blue")
plot(2*Re(z2),type="l",main="kzft(101,2)",xlab="t", ylab="y", ylim=c(-6,6))
lines(s,col="blue")
plot(2*Re(z3),type="l",main="kzft(101,3)",xlab="t", ylab="y", ylim=c(-6,6))
lines(s,col="blue")
par(mfrow=c(1,1))
}
