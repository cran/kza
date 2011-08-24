#I have a sine wave about period 10 with noise added to this signal N(0,16).  
#I randomly omitted 40% of the data to simulate a time series with missing data.  
#I denote percent missing by p.  Hence, p=.40.  

t <- 1:17000
x <- sin(2*pi*(1/10)*t)
#+rnorm(17000,0,4)
rand_idx <- sample(t,3400,replace=F)

x[rand_idx]<-NA
t[rand_idx]<-NA

tn<-na.omit(t)
r<-na.omit(x)

s<-1:17000
y<-sin(2*pi*(1/10)*s)
z.10<-kzft(r,m=200,k=3,f=.10,dim=1,index=tn)

z.sum<-2*Re(z.10)
plot(x[1:100],type='p',ylab="")
lines(z.sum[2:101],type='l',xlab="",ylab="")
lines(y[1:100],type="l",xlab="",ylab="",lty=3,lwd=3,col="red")
title(main="sin(2pi*(.1)t)+N(0,16) with reconstructed signal",sub="m=200,k=3,p=.2")
cor(z.sum[2:101],y[1:100])


2.) I then want to make the signal more complicated by considering a periodic signal that is the sum of two sine waves about periods 10 and 12.5.  I added noise N(0.16).  Then again randomly omitted data to consider a time series where p=.40.  I want to use the KZFT to reconstruct the signal.  I know that I have to apply KZFT twice about each of the known frequencies (.10 and .08)  I use kzft with m=200, k=3. 

t <- 1:17000
x <- sin(2*pi*(1/10)*t)+sin(2*pi*(1/12.5)*t)+rnorm(17000,0,4)
rand_idx <- sample(t,6800,replace=F)

x[rand_idx]<-NA
t[rand_idx]<-NA

tn<-na.omit(t)
r<-na.omit(x)

s<-1:17000
y<-sin(2*pi*(1/10)*s)+sin(2*pi*(1/12.5)*s)
z.10<-kzft(r,m=200,k=3,f=.10,dim=1,index=tn)
z.125<-kzft(r,m=200,k=3,f=.08,dim=1,index=tn)
z.sum<-2*Re(z.10)+2*Re(z.125)
plot(x[1:100],type='p',ylab="")
lines(z.sum[2:101],type='l',xlab="",ylab="")
lines(y[1:100],type='l',xlab="",ylab="",lty=3,lwd=3,col="black")
title(main="sin(2pi*(.1)t) +sin(2pi(.08)t)+N(0,16) with reconstructed signal",sub="m=200,k=3,p=.4")
cor(z.sum[2:101],y[1:100])

