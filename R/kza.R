kz <- function(x,q,k = 3) {
    if (length(dim(x)) > 3) stop("Too many dimensions.")
    storage.mode(x) <- "double"
    x <- .Call("kz", x, as.integer(round(q/2)), as.integer(k), PACKAGE="kza")
}

kza <- function(x,q,y=NULL,k=3,m=round(0.05*q),tol=1.0e-5) {
    if (length(dim(x)) > 3) stop("Too many dimensions.")
    if (is.null(y)) { 
    	y<-kz(x,q=q,k=k)
    }
    storage.mode(x) <- storage.mode(y) <- "double"
    x <- .Call("kza", x, y, as.integer(q), as.integer(k), as.integer(m), as.double(tol), PACKAGE="kza")
}

kzsv <- function(y,kza=NULL,kz=NULL,q,k=3,m=round(0.05*q),tol=1.0e-5) {
    w=round(q/2)
    if (is.null(kza) || is.null(kz)) {
		kz <- kz(y,q=w,k=k)
    	kza <- kza(y,q=w,k=k,m=m)
    }
    storage.mode(kza) <- storage.mode(kz) <- "double"
    s <- .Call("kzsv", kza, kz, as.integer(w), as.integer(k), as.integer(m), as.double(tol), PACKAGE="kza")
    
	par(mfrow=c(3,1))
	plot(y,type="l")
	plot(kza,type="l")
    plot(sqrt(s/mean(s))/1.96,ylab="sigma",type="l")
}

kzs <- function(y,m=round(length(y)/10),k = 1,t=NULL) {
	if (is.null(t)) a<-kzft(y,m=m,f=0,k=k,dim=1,trim=FALSE)[1:length(y)]
	else a<-kzft(y,m=m,t=t,f=0,k=k,dim=1)[1:length(y)]
	
	if (is.ts(y)) ans<-ts(2*Re(a),start=start(y),frequency=frequency(y))
	else ans<-2*Re(a)
	return (ans)
}
