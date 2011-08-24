#ckzft

#f<-0.03
#noise<- rnorm(length(t))

library(kza)

period=101
f<-1/period
t<-1:2000
s<-1*sin(2*pi*f*t)
x<-s

#a<-kzft(x,m=101,f=f,index=t)

rand_idx <- sample(t,1000,replace=F)

x[rand_idx]<-NA
t[rand_idx]<-NA

#x<-as.vector(na.omit(x))
#t<-as.vector(na.omit(t))

a<-kzft(x,m=51,f=f)



y.r<-vector(mode="numeric", length=length(t))
y.i<-vector(mode="numeric", length=length(t))
m=101
scale=1

xo=x
to=t

x.na<-ifelse (!!(t%%2)==TRUE, x, NA)
t.na<-ifelse (!!(t%%2)==TRUE, t, NA)

t<-as.vector(na.omit(t.na))
x<-as.vector(na.omit(x.na))

index<-t
s<-.C("ckzft", as.double(y.r), as.double(y.i), as.double(x),as.integer(length(x)), as.double(index), as.double(m), as.double(scale),as.double(f), NAOK=TRUE, DUP=FALSE, PACKAGE="kza")


#s<-.C("ckzft", as.double(y.r), as.double(y.i), as.double(x),as.integer(length(x)), as.double(t), as.double(m), as.double(scale),as.double(f), NAOK=TRUE, DUP=FALSE, PACKAGE="kza")




period=100
f<-1/period
t<-1:2000
s<-3*sin(2*pi*f*t)
x<-s

y.r<-vector(mode="numeric", length=length(t))
y.i<-vector(mode="numeric", length=length(t))
m=101
scale=1

xo=x
to=t

rand_idx <- sample(t,1000,replace=F)

x[rand_idx]<-NA
t[rand_idx]<-NA

x.na=x
t.na=t

#t[5]=NA
#x[5]=NA
#t[12]=NA
#x[12]=NA
t<-as.vector(na.omit(t))
x<-as.vector(na.omit(x))

s<-.C("ckzft", as.double(y.r), as.double(y.i), as.double(x), as.integer(length(x)), as.double(t), as.double(m), as.double(scale), as.double(f), NAOK=TRUE, DUP=FALSE, PACKAGE="kza")
plot(to,x.na,type='l',ylim=c(-4,4))
lines(2*y.r,col="blue")



m=101
f<-1/period
f<-f/2

s<-.C("ckzft", as.double(y.r), as.double(y.i), as.double(x), as.integer(length(x)), as.double(t), as.double(m), as.double(scale), as.double(f), NAOK=TRUE, DUP=FALSE, PACKAGE="kza")

plot(x,type='l')
lines(2*y.r,col="blue")
