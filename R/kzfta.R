.packageName <- "kza"

smooth<-function(p,c=0.01,w=length(p),method="DZ") 
{
    if (is.vector(p)) {
        if (is.null(w)) w<-length(p)
        func<-"smooth_kzp"
        N<-length(p)
        spg<-rep(0,N)
    } else { #assume kztp
        if (method != "DZ") 
            stop("Only method DZ is available for smoothing kztp at this time.")
        func<-"smooth_kztp"
        if (is.complex(p) == TRUE) {
            p<-Mod(p)
        }
        p<-as.matrix(p)
        if (is.null(w)) w<-dim(p)[1]
        N<-dim(p)[1]
        spg<-array(0,dim=c(N,N))
    }

    s <- .Call(func,
        p,
        as.double(c),
        as.integer(w),
        as.double(spg),
        as.character(method),
        PACKAGE="kza")
        
    if (is.matrix(p)) { dim(s)<-c(N,N) }
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

kzft <- function(x, m=100, k=1, f=NULL, dim=1)
{
    if (((m-1)*k+1)>length(x))
   	    stop("invalid 'm' or 'k':(m-1)*k+1 should be less equal length of the input vector x")

	k<-as.integer(k)
	m<-as.integer(m)
	if (is.null(nrow(x))) n<-length(x)-m+1 else n=nrow(x)-m+1
	z<-matrix(nrow=n, ncol=m, byrow=TRUE)
	for(i in 1:n) {
		z[i,]=rev(fft(x[i:(m+i-1)], inverse=TRUE))/m
	}

	if (is.null(f)) {
		s<-which.max(colMeans(abs(Re(z)))[1:(m/2)])
		f<-s/m
	} else {
		s<-f*m
	}
	k<-k-1
	if (k>0) { 
		z<-as.vector(z[,s])
		z<-kzft(z,m=m,k=k,f=f,dim=dim)$Complex
	}
	if (dim==1 & k==0) z<-as.vector(z[,s])
	lst<-list(Complex=z, freq=f)
    return(lst);
}

ff_kzft <- function(xr, xi=NULL, m=100, k=1, f=NULL, dim=1)
{
	if (is.null(xi)) { xi<-ff(vmode="double",length=length(xr))	}
	
	s<-f 
	k<-as.integer(k)
	m<-as.integer(m)
	n=nrow(xr)
	if (is.null(n)) {n=length(xr)-m+1} else n=n-m+1
	zr<-ff(vmode="double", dim=c(n,m))
	zi<-ff(vmode="double", dim=c(n,m))
	for(i in 1:n) {
		x<-complex(real=xr[i:(m+i-1)], imaginary=xi[i:(m+i-1)])
		z<-rev(fft(x, inverse=TRUE))
		zr[i,]=Re(z)/m
		zi[i,]=Im(z)/m
	}
	if (is.null(f)) {
		s<-which.max(colMeans(abs(zr[,])))
		f<-s/m
	} else {
		s<-f*m
	}

	k<-k-1
	if (k>0) {
		zr<-as.ff(zr[,s])
		zi<-as.ff(zi[,s])
		z<-ff_kzft(zr,zi,m,k,f=f,dim=dim)
		zr<-z$Re; zi<-z$Im;
	}
	if (dim==1 & k==0) {
		zr<-as.ff(zr[,s])
		zi<-as.ff(zi[,s])
	}
	z<-list(Re=zr, Im=zi, freq=f)
    return(z);
}

kzp <- function(z, m, k)
{
	if (!is.matrix(z)) stop("KZP requires matrix as input.")
	M<-(m-1)*k+1
	d<-apply(z,2,function(z) {(abs(z)^2)*M})
	if (is.null(dim(d))) a=d  else a<-colMeans(d)
    a<-a[1:round(m/2)]
    return(a)
}

kztp <- function(x, m, k, box=c(0,0.5,0,0.5))
{
    x<-as.vector(x)
    if (((m-1)*k+1)>length(x)) stop("The value of (m-1)*k+1 needs to be less then the length of the input data x")

	z<-kzft(x, m, k, dim=2)$Complex
	d<-kztp_conj(z, box)
    return(d)
}

kztp_conj <- function(z, box=c(0,0.5,0,0.5)) 
{
	if (!is.matrix(z)) stop("input must be a matrix.")
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
