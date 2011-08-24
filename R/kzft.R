.packageName <- "kza"

coeff<-function(m, k){as.vector(polynomial(rep(1,m))^k/m^k)}

kzft<-function (x, m = length(x), k = 1, f = NULL, dim = 1, index = NULL, trim=FALSE) 
{
	if (dim > 2) stop("kzft only supports up to 2 dimensions.")
	k <- as.integer(k)
	m <- as.integer(m)
	loops=k

	if (dim==1 & !is.null(f)) {
		
		scale=1
		if (is.null(index)) { 
			index=seq(1:length(x));
			index<-index[!is.na(x)];
		}

		if (max(is.na(index))) stop("Index cannot have NA values");
		for (i in 1:k) {
			x<-as.vector(na.omit(x))
			index<-as.vector(na.omit(index))
			y.r<-vector(mode="numeric", length=(max(index)));
			y.i<-vector(mode="numeric", length=(max(index)));
			s<-.C("ckzft", as.double(y.r), as.double(y.i), as.double(x), as.integer(length(x)), 
				as.double(index), as.double(m), as.double(scale), as.double(f), NAOK=TRUE, DUP=FALSE, PACKAGE="kza");
			x<-2*y.r;
			index<-seq(1:length(x));
			index<-index[!is.na(x)];
		}
		if (trim) { y.r<-y.r[(k*(m-1)/2):(length(y.r)-k*(m-1)/2)]; y.i<-y.i[(k*(m-1)/2):(length(y.i)-k*(m-1)/2)]; }
		z<-complex(real=y.r, imag=y.i);
	
		return (z)
	}

	    
	######
	# if no indexing set t and a frequency (f) is not supplied, assume a spectrum analysis
	#####
	if (is.null(index)) {
		if (is.null(f)) { 
			index<-seq(1,length(x)) 
		} else {
			if (f==0) { index<-seq(m/2,length(x)+m/2-1) }
			else { index<-seq(1,length(x)) }
		}
	}
	
	n<-max(index)
	z <- matrix(nrow = n, ncol = m, byrow = TRUE)
	z<-.Call("kzftwz",x,as.integer(index),as.integer(m),as.matrix(z))
	
	if (is.null(f)) {
		s <- which.max(colMeans(abs(Re(z)))[1:(m/2)])
		f <- (s-1)/m
	} else { s <- f * m + 1	}
	
	if (k>1) {
		if (f==0) { index<-seq(m/2,length(z[,s])+m/2-1) }
		else { index<-seq(1,length(z[,s])) }
		
		for (i in 2:k) {
			z<-.Call("kzftwz",z[,s],as.integer(index),as.integer(m),as.matrix(z))
		}
	}

	## remove m at the end of array
	if (trim & (n-m*loops)>0) z<-z[1:(n-m*loops),]
	if (dim==1) z<-as.vector(z[,s])
	    
	return (z)
}

Rkzft <- function(x, m=length(x), k=1, f=NULL, dim=2, c=NULL)
{
	k<-as.integer(k)
	m<-as.integer(m)

	if (is.null(c)) c<-coeff(m,k)
	if (is.null(nrow(x))) n<-length(x)-(k*(m-1)+1) else n=nrow(x)-m+1
	n<-max(1,n)
	z<-matrix(nrow=n, ncol=(k*(m-1)+1), byrow=TRUE)
	for(i in 1:n) {
		y<-x[i:(k*(m-1)+i)]*c
		z[i,]=fft(y)
	}

	if (is.null(f)) {
		s<-which.max(colMeans(abs(Re(z)), na.rm = TRUE)[1:(m/2)])
		f<-(s-1)/(k*(m-1)+1)
	} else {s<-round(f*(k*(m-1)+1)+1)}
	
	#if ((n-m*k)>0) z<-z[1:(n-m*k),]
	if (dim==1) z<-as.vector(z[,s])
    return(z);
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

kzp <- function(y, m=length(y), k=3, index=NULL, double_frequency=FALSE)
{
	M<-(m-1)*k+1
	z<-kzft(y,m=m,k=k,index=index,dim=2)
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

Rkzp <- function(y, m=length(y), k=3, f=NULL, double_frequency=FALSE)
{
	M<-(m-1)*k+1
	z<-Rkzft(y,m=m,k=k,f=f,dim=2)
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

periodogram <- function(y) {
	fourier<-fft(y)
	
	magnitude<-Mod(fourier)
	 
	# extract the phase which is atan(Im(fourier)/Re(fourier))
	phase<-Arg(fourier)
	
	# select only first half of vectors
	magnitude_firsthalf <- magnitude[1:(length(magnitude)/2)]
	phase_firsthalf<-phase[1:(length(magnitude)/2)]
	 
	# generate x-axis
	x.axis <- 1:length(magnitude_firsthalf)/length(magnitude)
	
	p<-cbind(x.axis, magnitude_firsthalf)
	return(p)
}

plot.kzp <- function(x, ...)
{
	if (is.null(x$smooth_periodogram)) dz<-x$periodogram else dz<-x$smooth_periodogram
	omega<-seq(0:(length(x$periodogram)-1))/(x$k*(x$window-1))
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
    print((mlist)/(object$k*(object$window)), digits=digits, ...)

    cat("\n Periods of interest:\n")
    print((object$k*(object$window))/(mlist), digits=digits, ...)
    invisible(object)
}

max_freq <- function(x, m, start=1, t=NULL)
{
    m <- as.integer(m)
    n <- length(x)
    
    if (is.null(t)) {t<-seq(1,length(x)) } 
    n<-max(t)

    z <- matrix(nrow = n, ncol = m, byrow = TRUE)
    z<-.Call("kzftwz",x,as.integer(t),as.integer(m),as.matrix(z))

	start=start
    s <- which.max(colMeans(abs(Re(z)))[start:(m/2)])
    f <- (s-1)/m
   	
    return (f)
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

nonlinearity.kzp<-function(x)
{
    n<-length(x)

    s<-rep(0,n)
    for (t in 2:(n-1)) {
        s[t]<-abs(x[t+1]-2*x[t]+x[t-1])
    }

    sq<-array(0, dim=c(n, n))

    for (i in (1:n)) for (j in (2:n)) {
        sq[i,j]<-sum(s[(max(1,(i-j+1))):(min(n,(i+j-1)))])
    }
    return(list(total=sum(s), matrix=sq))
}

variation.kzp<-function(x)
{
    n<-length(x)
    s=c(diff(x)^2,0)

    q<-array(0, dim=c(n, n))
    for (i in (1:n)) for (j in (2:n)) {
        q[i,j]<-sum(s[(max(1,(i-j+1))):(min(n,(i+j-2)))])
    }
    return(list(total=sum(s), matrix=q))
}

smooth.kzp<-function(object, log=TRUE, smooth_level=0.05, method = "DZ")
{
	if (class(object)!='kzp') stop ("Object type needs to be kzp.")
    n<-length(object$periodogram)
    spg<-rep(0,n)
    m<-rep(0,n)

	if (log==TRUE) p=log(object$periodogram) else p=object$periodogram
    if (method == "DZ") q<-variation.kzp(p)
    else if (method == "NZ") q<-nonlinearity.kzp(p)

    cc<-smooth_level*q$total

    for ( i in (1:n) ) {
        m[i]<-sum(q$matrix[i,1:n]<=cc)
        spg[i]<-mean(p[(max(1,(i-m[i]+1))):(min(n,(i+m[i]-1)))])
    }
    
    object$smooth_periodogram<-spg
    object$smooth_method=method
    return(object)
}

