coeff<-function(m, k){
      poly<-polynomial(rep(1,m))
      polyk<-poly^k
      coeff<-as.vector(polyk)
      coeff<-coeff/m^k
      return(coeff)
}

Rkzft <- function(x, m=length(x), k=1, f=NULL, dim=1, coeff=NULL)
{
	k<-as.integer(k)
	m<-as.integer(m)

	if (is.null(coeff)) coeff<-coeff(m,k)
	if (is.null(nrow(x))) n<-length(x)-(k*(m-1)+1) else n=nrow(x)-m+1
	z<-matrix(nrow=n, ncol=(k*(m-1)+1), byrow=TRUE)
	for(i in 1:n) {
		y<-x[i:(k*(m-1)+i)]*coeff
		z[i,]=fft(y)
	}

	if (is.null(f)) {
		s<-which.max(colMeans(abs(Re(z)))[1:(m/2)])
		f<-(s-1)/m
	} else {s<-f*m+1}
	
	if (!is.vector(z)) {
		if ((n-1)>0) z<-z[1:(n-1),]
		if (dim==1) z<-as.vector(z[,s])
	}
    
    return(z);
}

Rkzp <- function(y, m=length(y), k=1, f=NULL, double_frequency=FALSE)
{
	M<-(m-1)*k+1
	z<-Rkzft(y,m=m,k=k,f=f,dim=2)
	d<-apply(z,2,function(z) {(abs(z)^2)*M})
	if (is.null(dim(d))) a=d  else a<-colMeans(d)

	if (!double_frequency) a<-a[1:(m/2)]

	structure(list(
		periodogram = a,
		window=m,
		var=var(y),
		call=match.call()
            ),
        class = "kzp")
}

	if (smooth) {periodogram=smooth_kzp(log(a), c=smooth_level )} else {periodogram=a}

nonlinearity.kzp<-function(pg,K=length(pg)){
    N<-length(pg)

    S<-rep(0,N)
    S[1]<-0
    S[N]<-0
    for ( t in 2:(N-1) ) {
        S[t]<-abs(pg[t+1]-2*pg[t]+pg[t-1])
    }

    total<-sum(S)

    sq<-array(0, dim=c(N, K))

    for ( i in (1:N) ) for ( j in (2:K) ) {
        sq[i,j]<-sum(S[(max(1,(i-j+1))):(min(N,(i+j-1)))])
    }

    lst<-list(total=total, sqmatrix=sq)
    return(lst)
}

variation.kzp<-function(pg,K=length(pg)){
    N<-length(pg)
    S=c(diff(pg)^2,0)

    total<-sum(S)

    sq<-array(0, dim=c(N, K))
    for ( i in (1:N) ) for ( j in (2:K) ) {
        sq[i,j]<-sum(S[(max(1,(i-j+1))):(min(N,(i+j-2)))])
    }

    lst<-list(total=total, sqmatrix=sq)
    return(lst)
}

smooth.kzp<-function(pg,c,K=length(pg),method = "DZ") {
    N<-length(pg)
    spg<-rep(0,N)
    m<-rep(0,N)

    if      (method == "DZ") sq<-variation.kzp(pg,K)
    else if (method == "NZ") sq<-nonlinearity.kzp(pg,K)

    cc<-c*sq$total

    for ( i in (1:N) ) {
        m[i]<-sum(sq$sqmatrix[i,1:K]<=cc)
        spg[i]<-mean(pg[(max(1,(i-m[i]+1))):(min(N,(i+m[i]-1)))])
    }

    lst<-list(periodogram=spg,number=m)
    return(lst)
}


plot.kzp <- function(x, ...)
{
	if (is.null(x$smooth_periodogram)) dz<-x$periodogram else dz<-x$smooth_periodogram
	if (length(x$periodogram)==x$window) {to=1} else {to=1/2}
	omega<-seq(from=0,to=to,length.out=length(x$periodogram))

	plot(omega, dz, type="l", xlab="Frequency", ylab="")
}

summary.kzp <- function(object, digits = getOption("digits"), top=NULL, ...)
{
	cat(" Call:\n ")
	dput(object$call, control=NULL)

	if (is.null(object$smooth_periodogram))
		d<-diff(diff(object$periodogram))
	else
		d<-diff(diff(object$smooth_periodogram))
		
	if (is.null(top)) mlist=which(d< -2)
	else {
		mlist<-rep(0,top)
		for (i in 1:top) {
			mlist[i]<-which.max(-d)
			d[which.max(-d)]=NA			
		}
	}

    cat("\n Frequencies of interest:\n")
    print(mlist/(object$k*(object$window-1)), digits=digits, ...)

    cat("\n Periods of interest:\n")
    print((object$k*(object$window-1))/mlist, digits=digits, ...)
    invisible(object)
}


t<-1:1000
f<-0.03
m<-100
s<-3*sin(2*pi*f*t)
noise<- rnorm(length(t))
x<-s

system.time(z1<-kzft(x,m=m,k=1, dim=1))
system.time(z2<-kzft(x,m=m,k=2, dim=1))
system.time(z3<-kzft(x,m=m,k=3, dim=1))
system.time(z4<-kzft(x,m=m,k=4, dim=1))

par(mfrow=c(2,2))
plot(x[1:1000],type="l",main="Original time series",xlab="t", ylab="y")
plot(2*Re(z1),type="l",main="kzft(m,1)",xlab="t", ylab="y", ylim=c(-4,4))
lines(s,col="red")
plot(2*Re(z2),type="l",main="kzft(m,2)",xlab="t", ylab="y", ylim=c(-4,4))
lines(s,col="red")
plot(2*Re(z3),type="l",main="kzft(m,3)",xlab="t", ylab="y", ylim=c(-4,4))
lines(s,col="red")
par(mfrow=c(1,1))


t<-1:8000
f1<-0.03
f2<-0.04
noise<-35*rnorm(length(t))
s<-3*sin(2*pi*f1*t)+3*sin(2*pi*f2*t)
system.time(b<-kzp(s+noise,m=1000,k=3))

c<-smooth.kzp(b, smooth_level=0.01)

par(mfrow=c(3,1))
plot(periodogram(s+noise),type='l', main='Raw Periodogram')
plot(b)
plot(c)
par(mfrow=c(1,1))

max(b$periodogram)/(3*var(s+noise))

a<-log(b$periodogram)
y1<-smooth.kzp(a,0.01, method="NZ")$periodogram
y2<-smooth.kzp(a,0.01, method="DZ")$periodogram


t<-1:10000
f<-0.03
s<-3*sin(2*pi*f*t)
noise<- rnorm(length(t))
x<-s

#system.time(z<-kzft(x,m=100,k=3, dim=1))

system.time(z1<-kzft(x,m=100,k=3, dim=1))
system.time(zr<-Rkzft(x,m=100,k=3, dim=1))

par(mfrow=c(3,1))
plot(x[1:1000],type="l",main="Original time series",xlab="t", ylab="y")
plot(2*Re(z1)[1:1000],type="l",main="kzft(1000,1)",xlab="t", ylab="y", ylim=c(-4,4))
lines(s,col="red")
plot(2*Re(zr)[1:1000],type="l",main="kzft(1000,3)",xlab="t", ylab="y", ylim=c(-4,4))
lines(s,col="red")
par(mfrow=c(1,1))

plot(2*Re(zc)[1:1000],type="l",main="kzft(1000,2)",xlab="t", ylab="y", ylim=c(-4,4))
lines(s,col="red")

