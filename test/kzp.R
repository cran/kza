Rkzp <- function(y, m=length(y), k=3, double_frequency=FALSE)
{
	M<-(m-1)*k+1
	n<-length(y)
    z <- matrix(nrow = n, ncol = (m-1), byrow = TRUE)
	for (i in 1:(m-1)) {
		z[,i]<-Rkzft(y,m=m,k=k,f=i/m,dim=1,trim=FALSE)
	}
		
	d<-apply(z,2,function(z) {(abs(z)^2)*M})
	if (is.null(dim(d))) a=d  else a<-colMeans(d)

	if (!double_frequency) a<-a[1:(m/2)]

	structure(list(
		periodogram = a,
		window=m,
		k=k,
		var=var(y),
		smooth_periodogram=NULL,
		smooth_method=NULL,
		call=match.call()
            ),
        class = "kzp")
}

x <- 1:4
fft(x)
fft(fft(x), inverse = TRUE)/length(x)

r<-c(1,rep(0,7))
fft(r)
r<-c(0,1,rep(0,6))
fft(r)


t<-1:7000
f1<-0.03
f2<-0.04
noise<-rnorm(length(t),0,5)
s<-3*sin(2*pi*f1*t)+3*sin(2*pi*f2*t)
a<-kzp(s,m=500,k=1)
par(mfrow=c(3,1))
plot(kzp(s+noise),main="KZP periodogram (smooth=FALSE)")
plot(kzp(s+noise,smooth=TRUE,smooth_level=0.05))
plot(periodogram(s+noise),type='l')
par(mfrow=c(1,1))

summary(a, top=2)


t<-1:6000
f1<-0.030
f2<-0.039
noise<-10*rnorm(length(t))
amp=3
s<-amp*sin(2*pi*f1*t)+amp*sin(2*pi*f2*t)

system.time(b<-kzp(s+noise,m=500,k=3))
c<-smooth.kzp(b, smooth_level=0.1)
a<-kzpi(s+noise,m=500,k=3)

par(mfrow=c(4,1))
plot(periodogram(s+noise),type='l', main='Raw Periodogram')
plot(a)
plot(b)
plot(c)
par(mfrow=c(1,1))


t<-1:6000
f1<-0.03
f2<-0.04
noise<-15*rnorm(length(t))
amp=1.5
s<-amp*sin(2*pi*f1*t)+amp*sin(2*pi*f2*t)
system.time(a<-kzp(s+noise,500,k=3))
b<-smooth.kzp(a, smooth_level=0.08)
par(mfrow=c(3,1))
plot(periodogram(s+noise),type='l')
plot(a)
plot(b)
par(mfrow=c(1,1))

# signal/noise
signal<-kzft(s+noise,m=500,k=3,dim=1)
print(sqrt(var(2*Re(signal))/var(s+noise)))
