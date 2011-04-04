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
### linear trend and 12 month cycle
###
yrs <- 20
t <- seq(0,yrs,length=yrs*365)
amp<-10
y <- amp*sin(2*pi*t*12+10)
y<-ts(y,start=c(1980,10),frequency=365)
signal<-y

m<-seq(1, amp, length.out=length(y))
y<-y+m
set.seed(6); e <- rnorm(n = length(t), sd = 5*amp)
SNR<-mean(abs(y)^2)/(mean((abs(e)^2)))
y<-y+e
a<-decompose.kz(y)
plot.decompose.kz(a)
#plot(a$one_season)

xstl<-stl(y,"per")
plot(xstl)

r1<-signal+m-a$fitted
r3<-signal+m-(xstl$time.series[,1]+xstl$time.series[,2])

var(r1)
var(na.omit(r3))

plot(cbind(
y,
signal+m,
a$fitted,
(xstl$time.series[,1]+xstl$time.series[,2]))
)

p<-predict(a,n.ahead=365)
plot(append(a$fitted,p)

###
### exp trend and monthly cycle
###
yrs <- 20
t <- seq(0,yrs,length=yrs*365)
amp<-10
y <- amp*sin(2*pi*t*12+10)
y<-ts(y,start=c(1980,10),frequency=365)
signal<-y

tr<-seq(1, amp/2, length.out=7300)
y.orig<-signal+exp(tr)+300
y<-y.orig
set.seed(6); e <- rnorm(n = length(t), sd = 50)
SNR<-mean(abs(y)^2)/(mean((abs(e)^2)))
y<-log(y+e)

a<-decompose.kz(y)
plot.decompose.kz(a)
a$seasonal_trend<-ts(exp(a$seasonal_trend)*100-100,start=start(a$trend),frequency=frequency(a$trend))
a$seasonal_mean<-ts(exp(a$seasonal_mean)*100-100,start=start(a$trend),frequency=frequency(a$trend))
a$one_season<-ts(exp(a$one_season)*100-100,start=start(a$trend),frequency=frequency(a$trend))
a$trend<-ts(exp(a$trend),start=start(a$trend),frequency=frequency(a$trend))


plot(a$trend)
lines(ts(exp(tr)+300,start=start(signal),frequency=365),col="red")

plot(signal)
lines(a$seasonal_trend,col="blue")

plot(exp(a$fitted))
lines(y.orig,col="blue")

plot.decompose.kz(a)

xd<-decompose(y)
plot(xd)

xstl<-stl(y,"per")
plot(xstl)

r1<-y.orig-exp(a$fitted)
r2<-y.orig-exp(xd$seasonal+xd$trend)
r3<-y.orig-exp(xstl$time.series[,1]+xstl$time.series[,2])

var(r1)
var(na.omit(r2))
var(na.omit(r3))

plot(y.orig)
lines(exp(a$fitted),col="blue")

plot(y.orig)
lines(exp(xd$seasonal+xd$trend),col="blue")

plot(y.orig)
lines(exp(xstl$time.series[,1]+xstl$time.series[,2]),col="blue")
