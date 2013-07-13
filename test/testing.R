#######
# this is an example of detection of a break point in a time series
#######
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
k.kza <- kza(y,q,y=k.kz,m=10)

par(mfrow=c(2,1))
plot(y,type="l", ylim=c(-3,3))
plot(signal,type="l",ylim=c(-3,3), main="Signal and KZA Reconstruction")
lines(k.kza$kza, col=4)

x <- c(rep(0,4000),rep(0.5,2000),rep(0,4000))
noise <- rnorm(n = 10000, sd = 1.0) # normally-distributed random variates
v <- x + noise
a<-kzsv(v, q=1000, k=3)

    
	par(mfrow=c(3,1))
	plot(y,type="l")
	plot(kza,type="l")
    plot(sqrt(s/mean(s))/1.96,ylab="sigma",type="l")


for(i in 1:10) {
       l<-l+years[i]
       m<- which.max(x[l:(l+years[(i+1)])]) - 1
       a<- s[l+m]
       p<- l+m
       peaks[i]<-p
       peaks.text[i]<-paste(substr(a,6,10))

       m<- which.min(x[l:(l+years[(i+1)])]) - 1
       a<- s[l+m]
       p<- l+m
       troughs[i]<-p
       troughs.text[i]<-paste(substr(a,6,10))
}

peaks <- function(x)
{
	# second derivative
	a<-diff(diff(x))	
	peaks=NULL
	i=0

	p<-a[a<(-1)]
	peaks<-rep(0,length(p))
	i=1
	for (j in p) {
		peaks[i] <-	which(a==j)
		i=i+1
	}
	
}

cluster <- function(x, span=15)
{
	p=NULL
	m=NULL
	i=1
	z=x[i]
	for (y in x) {
		if (abs(y-z)<span) p=c(p,y) else  {m=c(m,round(mean(p))); p=c(y); } 
		z=x[i]
		i=i+1
	}
	m=c(m,round(mean(p)))
}


peaks <- function(x, sigma=3, span=15)
{
	a<-x-mean(x)
	p2=NULL
	i=0

	sigma=sigma*sqrt(var(x))
	p<-a[a>sigma]
	p2<-rep(0,length(p))
	i=1
	for (j in p) {
		p2[i] <-	which(a==j)
		i=i+1
	}
	return (cluster(p2, span))
}


y <- c(rep(0,4000),rep(0.5,2000),rep(0,4000))
noise <- rnorm(n = 10000, sd = 1.0) # normally-distributed random variates
v <- y + noise
z<-kzsv(v, q=1000, k=3)

x<-sqrt(z$kzsv/mean(z$kzsv))
