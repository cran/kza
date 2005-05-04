kz <- function(x,q,k = 3) {
    if (length(dim(x)) > 2) stop("Too many dimensions.")
    x <- .Call("kz", x, q, k, NAOK=TRUE, PACKAGE="kza")
}

kza <- function(x,q,k=3,m=0,tol=1.0e-5) {
    if (length(dim(x)) > 2) stop("Too many dimensions.")
    z <- .Call("kz", x, q, k, NAOK=TRUE, PACKAGE="kza")
    x <- .Call("kza", x, z, q, k, m, tol, NAOK=TRUE, PACKAGE="kza")
}

kzsv <- function(x,q,k=3,m=0,tol=1.0e-5) {
    z <- .Call("kz", x, q, k, NAOK=TRUE, PACKAGE="kza")
    k <- .Call("kzsv", x, z, q, k, m, tol, NAOK=TRUE, PACKAGE="kza")
}

