t <- 1:17000
x <- sin(2*pi*(1/12)*t)

 z.12<-kzft(x,m=621,k=12,f=414/4968,dim=1,index=t, alg="C")
 z.sum<-2*Re(z.12)

 split.screen(c(1,1))
 plot(7105:7249,x[7105:7249],lty=4,type="b",col="red",xlab="time",ylab="")
 screen(1,FALSE)
 plot(7105:7249,z.sum[7105:7249],col="blue",type="b",xlab="",ylab="",axes=FALSE)
 axis(side=4)
 close.screen(all=TRUE)

t <- 1:3000
t2<- 1:3000
x <- sin(2*pi*(1/12)*t)
m=621
f=1/12

z.12<-kzft(x,m=m,k=8,f=f,dim=1,index=t2,alg="C")
z.sum<-2*Re(z.12)

split.screen(c(1,1))
plot(1000:1100,x[1000:1100],lty=4,type="b",col="red",xlab="time",ylab="")
screen(1,FALSE)
plot(1000:1100,z.sum[1000:1100],col="blue",type="b",xlab="",ylab="",axes=FALSE)
axis(side=4)
close.screen(all=TRUE)


t <- 1:3000
t2<- 1:3000
f=1/11
#x <- sin(2*pi*f*t)
x<-cos(2*pi*f*t)
m=11
k=3

z.12<-kzft(x,m=m,k=k,f=f,dim=1,index=t2,alg="C")
z.sum<-2*Re(z.12)

z.sum<-x
for (i in 1:3) {
	z.12<-kzft(z.sum,m=m,k=1,f=f,dim=1,index=t2,alg="C")
	z.sum<-2*Re(z.12)
}


split.screen(c(1,1))
plot(1000:1100,x[1000:1100],lty=4,type="b",col="red",xlab="time",ylab="")
screen(1,FALSE)
plot(1000:1100,z.sum[1000:1100],col="blue",type="b",xlab="",ylab="",axes=FALSE)
axis(side=4)
close.screen(all=TRUE)

a[1]=x[100-5]*exp(complex(r=0,i=-2*pi*f*(100-5)))
a[1]=a[1]+x[100-4]*exp(complex(r=0,i=-2*pi*f*(100-4)))
a[1]=a[1]+x[100-3]*exp(complex(r=0,i=-2*pi*f*(100-3)))
a[1]=a[1]+x[100-2]*exp(complex(r=0,i=-2*pi*f*(100-2)))
a[1]=a[1]+x[100-1]*exp(complex(r=0,i=-2*pi*f*(100-1)))
a[1]=a[1]+x[100-0]*exp(complex(r=0,i=-2*pi*f*(100-0)))
a[1]=a[1]+x[100+1]*exp(complex(r=0,i=-2*pi*f*(100+1)))
a[1]=a[1]+x[100+2]*exp(complex(r=0,i=-2*pi*f*(100+2)))
a[1]=a[1]+x[100+3]*exp(complex(r=0,i=-2*pi*f*(100+3)))
a[1]=a[1]+x[100+4]*exp(complex(r=0,i=-2*pi*f*(100+4)))
a[1]=a[1]+x[100+5]*exp(complex(r=0,i=-2*pi*f*(100+5)))

a<-rep(0,3000)
for (t in 100:200) {
	a[t]=0;
	for(j in -5:5) {
		a[t]=a[t]+x[t+j]*exp(0-j*pi*f*2i)
	}
	a[t]=a[t]/m;
}

f=1/11
m=11
a<-rep(0,3000)
a[100]=(1/m)*cos(2*pi*f*(100-5))*exp(0-2*pi*f*(-5i))
a[100]=a[100]+(1/m)*cos(2*pi*f*(100-4))*exp(0-2*pi*f*(-4i))
a[100]=a[100]+(1/m)*cos(2*pi*f*(100-3))*exp(0-2*pi*f*(-3i))
a[100]=a[100]+(1/m)*cos(2*pi*f*(100-2))*exp(0-2*pi*f*(-2i))
a[100]=a[100]+(1/m)*cos(2*pi*f*(100-1))*exp(0-2*pi*f*(-1i))
a[100]=a[100]+(1/m)*cos(2*pi*f*(100))*exp(0-2*pi*f*(0i))
a[100]=a[100]+(1/m)*cos(2*pi*f*(100+1))*exp(0-2*pi*f*(1i))
a[100]=a[100]+(1/m)*cos(2*pi*f*(100+2))*exp(0-2*pi*f*(2i))
a[100]=a[100]+(1/m)*cos(2*pi*f*(100+3))*exp(0-2*pi*f*(3i))
a[100]=a[100]+(1/m)*cos(2*pi*f*(100+4))*exp(0-2*pi*f*(4i))
a[100]=a[100]+(1/m)*cos(2*pi*f*(100+5))*exp(0-2*pi*f*(5i))


	



a[1]=a[1]+x[3]*exp(complex(r=0,i=-2*pi*f*2))
a[1]=a[1]+x[4]*exp(complex(r=0,i=-2*pi*f*3))
a[1]=a[1]+x[5]*exp(complex(r=0,i=-2*pi*f*4))
a[1]=a[1]+x[6]*exp(complex(r=0,i=-2*pi*f*5))






a[2]=x[2]*exp(r=0,i=

a[2]=exp(complex(r=0,i=4*pi/6))*x[3]
a[3]=exp(complex(r=0,i=3*pi/6))*x[4]
a[4]=exp(complex(r=0,i=2*pi/6))*x[5]
a[5]=exp(complex(r=0,i=1*pi/6))*x[6]
a[6]=exp(complex(r=0,i=0*pi/6))*x[7]
a[7]=exp(complex(r=0,i=-1*pi/6))*x[8]
a[8]=exp(complex(r=0,i=-2*pi/6))*x[9]
a[9]=exp(complex(r=0,i=-3*pi/6))*x[10]
a[10]=exp(complex(r=0,i=-4*pi/6))*x[11]
a[11]=exp(complex(r=0,i=-5*pi/6))*x[12]


## notes
z.12<-kzft(x,m=m,k=k,f=f,dim=1,index=t2,alg="C")
This has a phase shift when f and m are odd, but does not when even.


