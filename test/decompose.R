#SNR<-mean(abs(signal)^2)/(mean((abs(et)^2)))

yrs <- 20
t <- seq(0,yrs,length=yrs*365)
#y <- sin(2*pi*t) + sin(3*pi*t)

y <- 1000*sin(2*pi*t*12+10)
y<-ts(y,start=c(1980,10),frequency=365)
k=3

y <- 1000*sin(2*pi*t*12+10)
y<-ts(y,start=c(1980,10),frequency=365)

set.seed(6); e <- rnorm(n = length(t), sd = 2.0)
y<-y+e
a<-decompose.kz(y)

###
### Trend and 12 month cycle
###
amp<-10
y <- amp*sin(2*pi*t*12+10)
y<-ts(y,start=c(1980,10),frequency=365)
signal<-y

m<-seq(1, amp, length.out=7300)
y<-y+m
set.seed(6); e <- rnorm(n = length(t), sd = 3*amp)
SNR<-mean(abs(y)^2)/(mean((abs(e)^2)))
y<-y+e
a<-decompose.kz(y)
plot.decompose.kz(a)
plot(a$one_season)

xd<-decompose(y)
plot(xd)

xstl<-stl(y,"per")
plot(xstl)

r1<-signal+m-a$fitted
r2<-signal+m-(xd$seasonal+xd$trend)
r3<-signal+m-(xstl$time.series[,1]+xstl$time.series[,2])

var(r1)
var(na.omit(r2))
var(na.omit(r3))


###
### exp trend and monthly cycle
###
amp<-10
y <- amp*sin(2*pi*t*12+10)
y<-ts(y,start=c(1980,10),frequency=365)
signal<-y

m<-seq(1, amp/2, length.out=7300)
y<-y+exp(m)
set.seed(6); e <- rnorm(n = length(t), sd = 50)
SNR<-mean(abs(y)^2)/(mean((abs(e)^2)))
y<-y+e
a<-decompose.kz(y)
plot.decompose.kz(a)
plot(a$one_season)

xd<-decompose(y)
plot(xd)

xstl<-stl(y,"per")
plot(xstl)

r1<-signal+exp(m)-a$fitted
r2<-signal+exp(m)-(xd$seasonal+xd$trend)
r3<-signal+exp(m)-(xstl$time.series[,1]+xstl$time.series[,2])

var(r1)
var(na.omit(r2))
var(na.omit(r3))

