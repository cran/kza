kz <- function(x,q,k = 3) {
    if (length(dim(x)) > 3) stop("Too many dimensions.")
    storage.mode(x) <- "double"
    x <- .Call("kz", x, q, k, NAOK=TRUE, PACKAGE="kza")
}

kza <- function(x,q,kz=NULL,k=3,m=0,tol=1.0e-5) {
    if (length(dim(x)) > 3) stop("Too many dimensions.")
    storage.mode(x) <- storage.mode(kz) <- "double"
    if (!length(kz)) { 
        print("Determine moving average with KZ.")
        kz <- .Call("kz", x, q, k, NAOK=TRUE, PACKAGE="kza") 
    }
    x <- .Call("kza", x, kz, q, k, m, tol, NAOK=TRUE, PACKAGE="kza")
}

kzsv <- function(y=NULL,kza=NULL,kz=NULL,q,k=3,m=0,tol=1.0e-5) {
    storage.mode(y) <- "double"
    storage.mode(kza) <- storage.mode(kz) <- "double"
    if (!length(y) && !length(kza)) { 
        print ("Either y or kza must be specified.")
        return( NULL )
    }
    if (!length(kz)) { 
        print ("finding moving average with kz")
        kz <- .Call("kz", y, q, k, NAOK=TRUE, PACKAGE="kza") 
        kza <- .Call("kza", y, kz, q, k, m, tol, NAOK=TRUE, PACKAGE="kza") 
    }
    s <- .Call("kzsv", kza, kz, q, k, m, tol, NAOK=TRUE, PACKAGE="kza")
    return( s/mean(s))
}
