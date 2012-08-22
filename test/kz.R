library(kza)

f=1/365.25
t<-1:(365.25*10)
s<-50*sin(2*pi*f*t)
x<-ts(s,start=c(2001,1),freq=365.25)+seq(1:length(t))

plot(x)

## 2 dimensional
a <- matrix(rep(0,100*100),nrow=100)
a[35:70,35:70]<-1
a <- a + matrix(rnorm(100*100,0,1),nrow=100)
m = c(20,5)
k=3

z <- .Call("kz", a, as.vector(as.integer(floor(m))), as.integer(k))
x <- seq(1,100)
y <- x
op <- par(bg = "white")
persp(x, y, z, zlab="z", ticktype="detailed", theta = 60, phi = 45, col = "lightblue")


#y<-kz(a,m=m,k=1)
v <- kza(a,m=m,y=z,k=3,impute_tails=TRUE)
x <- seq(1,100)
y <- x
op <- par(bg = "white")
persp(x, y, v$kza, zlab="z", zlim=c(-0.5,2), ticktype="detailed", theta = 30, phi = 30, col = "lightblue")

par(mfrow=c(1,2))
image(a,col=gray(seq(0,1,1/255)))
image(v$kza,col=gray(seq(0,1,1/255)))
par(mfrow=c(1,1))

## 3 dimensional
w<-array(data=0, dim=c(100,100,100))
for (i in 35:65) {
    w[i:65,35:65,i] <- 1
}
w<-w+rnorm(n = 100*100*100, sd = 2)
m=c(20,10,5)
system.time(a<-kz(w,m,3))

x11()
for(i in 1:50) {
    image(matrix(b$kz[,,i],50,50),col=gray(seq(0,min(max(b$kz[,,i],na.rm=TRUE),1),1/255)))
}

system.time(b<-kza(w,y=a,m=m,k=3,impute_tails=TRUE))

######################
# movie of filtered object
#####################
x11()
for(i in 1:50) {
    image(matrix(b$kza[,,i],50,50),col=gray(seq(0,min(max(b$kza[,,i],na.rm=TRUE),1),1/255)))
}
