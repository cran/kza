
kza <- function(v,q,k = 3)
 .C("kza",
 as.double(v),
 as.integer(length(v)),
 as.integer(q),
 as.integer(k), 
 NAOK=TRUE, 
 PACKAGE="kza")

kz <- function(v,q,k = 3)
 .C("kz",
 as.double(v),
 as.integer(length(v)),
 as.integer(q),
 as.integer(k), 
 NAOK=TRUE,
 PACKAGE="kza")


kzsv <- function(v,q,f)
 .C("kzsv",
 as.double(v),
 as.integer(length(v)),
 as.integer(q),
 as.double(f), 
 NAOK=TRUE,
 PACKAGE="kza")
 
 
kzf <- function(v,q,d,k = 3)
 .C("kzf",
 as.double(v),
 as.integer(length(v)),
 as.integer(q),
 as.double(d),
 as.integer(k),
 PACKAGE="kza")
 
