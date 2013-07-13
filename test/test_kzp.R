N=1024
s<-rep(5,N)
fft(s)

s<-rep(0,N)
s[1]=1
fft(s)

x<-seq(0,1023,by=0.5)
y<-2*sin(2*pi*.5*x) #amplitude =2, frequency=0.5

t<-1:1000
f<-10/1000
x<-cos(2*pi*f*t)

# first harmonic
t<-1:1024
f<-1/1024
s<-cos(2*pi*f*t)
plot(t,s,type='l')
periodogram(s)
