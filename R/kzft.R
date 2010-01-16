.packageName <- "kza"

Rkzft <- function(x, m=100, k=1, f=NULL, dim=1)
{
	k<-as.integer(k)
	m<-as.integer(m)

	if (is.null(nrow(x))) n<-length(x)-m+1 else n=nrow(x)-m+1
	z<-matrix(nrow=n, ncol=m, byrow=TRUE)
	for(i in 1:n) {
		z[i,]=fft(x[i:(m+i-1)])/m
	}

	if (is.null(f)) {
		s<-which.max(colMeans(abs(Re(z)))[1:(m/2)])
		f<-s/m
	} else {
		s<-f*m+1
	}
	k<-k-1
	if (k>0) { 
		z<-as.vector(z[,s])
		z<-Rkzft(x=z,m=m,k=k,f=f,dim=dim)
	}
	if (dim==1 & k==0) z<-as.vector(z[,s])
    return(z);
}

##
# kzft
# y is data
# m is filter window
# t is indexing set
# k is the number of iterations
# dim is how many dimensions to return 
# 	currently the maximum number of dimensions is limited to two
##
kzft<-function (x, m, k = 1, f = NULL, dim = 2, t = NULL, trim=TRUE) 
{
	if(.Call("check_fftw")) {
	
		if (dim > 2) stop("kzft only supports up to 2 dimensions.")
		k <- as.integer(k)
	    m <- as.integer(m)
	    n <- length(x)
	    
		######
		# if no indexing set t and a frequency (f) is not supplied, assume a spectrum analysis
		#####
	    if (is.null(t)) {
			if (is.null(f)) { 
				t<-seq(1,length(x)) 
			} else {
				if (f==0) { t<-seq(m/2,length(x)+m/2-1) }
				else { t<-seq(1,length(x)) }
			}
		}
	
		n<-max(t)
	    z <- matrix(nrow = n, ncol = m, byrow = TRUE)
	    z<-.Call("kzftwz",x,as.integer(t),as.integer(m),as.matrix(z))
	
		if (is.null(f)) {
			s <- which.max(colMeans(abs(Re(z)))[1:(m/2)])
	    	f <- (s-1)/m
		} else { s <- f * m + 1	}
	
		if (k>1) {
			if (f==0) { t<-seq(m/2,length(z[,s])+m/2-1) }
			else { t<-seq(1,length(z[,s])) }
			
		    for (i in 2:k) {
				z<-.Call("kzftwz",z[,s],as.integer(t),as.integer(m),as.matrix(z))
			}
		}

		## remove m at the end of array
		if (trim) z<-z[1:(n-m),]
	
	    if (dim==1) z<-z[,s]
	} else {
		print("FFTW3 is highly recommend, see README.")
		z <- Rkzft(x=x, m=m, k = k, f = f, dim = dim)
	}
	
    return (z)
}

Rkzp <- function(y, m=round(length(y)/10), k=1, f=NULL)
{
	M<-(m-1)*k+1
	z<-Rkzft(y,m=m,k=k,f=f,dim=2)
	d<-apply(z,2,function(z) {(abs(z)^2)*M})
	if (is.null(dim(d))) a=d  else a<-colMeans(d)
#    a<-a[1:round(m/2)]
    return(a)
}

kzp <- function(y, m=round(length(y)/10), k=1, f=NULL, double_frequency=FALSE)
{
	if(.Call("check_fftw")) {
		z<-vector(mode="numeric", length=m)
	
		if (k==1) z<-.Call("kzp_1k",as.double(y),as.integer(m),as.double(z))
		else {
			if (is.null(f)) {
				z<-.Call("kzp_1k",as.double(y),as.integer(m),as.double(z))
				s <- which.max(abs(Re(z))[1:(m/2)])
			    	f <- (s-1)/m
			} else { 
				s <- f * m + 1
			}
			y<-kzft(x=y,m=m,k=(k-1),f=f,dim=1)
			z<-.Call("kzp",as.complex(y),as.integer(m),as.double(z))
		}
	} else { z<-Rkzp(y,m=m,k=k,f=f) }

	if (!double_frequency) z<-z[1:(m/2)]
	
	return (log(z))
}

kztp <- function(x, m, k, box=c(0,0.5,0,0.5))
{
    x<-as.vector(x)
    if (((m-1)*k+1)>length(x)) stop("The value of (m-1)*k+1 needs to be less then the length of the input data x")

	z<-kzft(x, m=m, k=k, dim=2)
	T<-dim(z)[1]
	m<-dim(z)[2]
	rp1=box[1]; rp2=box[2]; cp1=box[3]; cp2=box[4];
	rm1<-max(round(m*rp1),1)
    rm2<-round(m*rp2)
    cm1<-max(round(m*cp1),1)
    cm2<-round(m*cp2)
    delta.rm<-rm2-rm1+1
    delta.cm<-cm2-cm1+1
	zp<-array(NA,dim=c(delta.rm,delta.cm,T))
    for ( i in (1:delta.rm) ) for ( j in (1:delta.cm) ){
    	 zp[i,j,]<-z[,i+rm1-1]*z[,j+cm1-1]*Conj(z[,i+j+rm1+cm1-2])*m^2
    }               
    d<-rowMeans(zp,dim=2)
    return(d)
}

frequency_kzft <- function(x, m, start=0, t=NULL)
{
    m <- as.integer(m)
    n <- length(x)
    
    ######
    # if no indexing set t and a frequency (f) is not supplied, assume a spectrum analysis
    #####
    if (is.null(t)) {t<-seq(1,length(x)) } 
    n<-max(t)

    if(.Call("check_fftw")) {
        z <- matrix(nrow = n, ncol = m, byrow = TRUE)
        z<-.Call("kzftwz",x,as.integer(t),as.integer(m),as.matrix(z))
    } else {
		z<-Rkzft(x,m=m,dim=2)
    }

	start=start+1
    s <- which.max(colMeans(abs(Re(z)))[start:(m/2)])
    f <- (s-1)/m
   	
    return (f)
}

smooth_kzp<-function(x, c=0.01, w=length(x), method="DZ") 
{
    if (is.vector(x)) {
        if (is.null(w)) w<-length(x)
        func<-"smooth_kzp"
        N<-length(x)
        spg<-rep(0,N)
    } else { #assume kztp
        if (method != "DZ") 
            stop("Only method DZ is available for smoothing kztp at this time.")
        func<-"smooth_kztp"
        if (is.complex(x) == TRUE) {
            p<-Mod(x)
        }
        x<-as.matrix(x)
        if (is.null(w)) w<-dim(x)[1]
        N<-dim(x)[1]
        spg<-array(0,dim=c(N,N))
    }

    s <- .Call(func,
        x,
        as.double(c),
        as.integer(w),
        as.double(spg),
        as.character(method),
        PACKAGE="kza")
        
    if (is.matrix(x)) { dim(s)<-c(N,N) }
    return(s)
}

transfer_function <- function(m, k, lamda=seq(-0.5,0.5,by=0.01), omega=0)
{
      lamda<-lamda*2*pi
      omega<-omega*2*pi

      N<-length(lamda)
      tf<-array(0,dim=c(N, m))

      for ( j in (1:m) ){
         tf[,j]<-exp(1i*(lamda-omega)*j)
      }

      tf1<-rep(0,N)
      for ( i in (1:N) ){
         tf1[i]<-sum(tf[i,])
      }

      tf2<-(1/m)*tf1
      tf2<-abs(tf2)^k
      return(tf2)
}

plot.kzp <- function(x, c=0.01, w=length(x), method="DZ", smooth_only=FALSE, ...)
{
	dz<-smooth_kzp(x, c, w, method="DZ")

	n=length(x)
	omega<-seq(0,1/2,length=n)

	if (smooth_only) plot(omega, dz, main="Smoothed Periodogram DZ method", type="l", xlab="Frequency", ylab="")
	else {
		par(mfrow=c(2,1))
		plot(omega, x, main="Raw periodogram", type="l", xlab="Frequency", ylab="")
		plot(omega, dz, main="Smoothed Periodogram DZ method", type="l", xlab="Frequency", ylab="")
	}
}
