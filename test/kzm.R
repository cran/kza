library(chron)

ny<-subset(fars, STATE==36)
ny<-subset(ny,!is.na(ny$latitude) | !is.na(ny$longitud))

d<-paste(ny$YEAR,ny$MONTH,ny$DAY,sep="/")
d<-chron(d,format="y/m/d")
d<-d-min(d, na.rm=TRUE)+1
f<-ny$FATALS

dx=100
dy=100
a<-cbind(cut(ny$longitud,dx), cut(ny$latitude,dy), d, ny$FATALS)
b<-data.frame(a)
names(b)<-c("V1","V2","V3","V4")
nams<-c("V1","V2","V3")
s = b[do.call(order, b[nams]), ]


df<-b
# subset
which((df[,1]>=1 & df[,1]<=3) & (df[,2]>=41 & df[,2]<=46) & (df[,3]>=1 & df[,3]<=365))
which((df[,1] %in% 1:3) & (df[,2] %in% 41:46) & (df[,3] %in% 1:365))
subset(df, V1 %in% 1:3 & V2 %in% 41:46 & V3 %in% 1:365)

dx=3
dy=3
dz=365

kzm<-function(df) {
	a<-array(NA, dim=c(100,100,3294))
	for (k in 1:100) {
		for (j in 1:(100-dy)) {
			for (i in 1:(100-dx)) {
				#a[i,j,k]=mean(df[which((df[,1] %in% i:(i+3)) & (df[,2] %in% j:(j+3)) & (df[,3] %in% k:(k+365)),4])
				a[i,j,k]=mean(subset(df, V1 %in% i:(i+dx) & V2 %in% j:(j+dy) & V3 %in% k:(k+dz))[4])
			}
		}
	}
	return (a)	
}
	
	
j=42	
for (i in 1:(100-dx)) {
	a[i,j,1]=mean(df[which((df[,1] %in% i:(i+3)) & (df[,2] %in% j:(j+3)) & (df[,3] %in% 1:365)),4])
}

kzm <- function(z, data, window)
{
	m=max(z);
	if (is.array(z)) {
		n=length(dim(z)) 
		d=dim(z)
	} else {
		n=1
		d=length(data)
	}
	out=vector(length=m^n,mode="numeric")
	wts=1
	window=c(3,3)
	
	.C("kz_ma", out,  z, length(z), data, d, n, window, wts, NAOK=TRUE, DUP=FALSE, PACKAGE="kza") 

//	s <- .Call("kzm", a, length(a), z, length(z), data, d, n, window, wts, PACKAGE="kza")	
}



i1<-c(1,1,1,2,2,2,3,3,4,5,5,5)
i2<-c(1,2,4,1,2,3,1,3,1,1,2,3)
im<-cbind(i1,i2)
data<-c(3,2,1,6,4,1,2,3,1,1,1,3)
window<-c(3,3)

x<-im
w<-c(1)

w=1
dx=dim(im)
n<-as.integer(length(dx))
x<-as.vector(im)
m=c(5,4);
len=5*4
a=vector(length=len,mode="numeric")

s<-.C("kz_ma", as.double(a), as.double(m), as.double(x), as.integer(length(x)), as.double(data), as.double(dx), as.integer(length(dx)), as.double(window), as.double(w), NAOK=TRUE, DUP=FALSE, PACKAGE="kza")






x<-c(1,1,1,2,2,2,3,3,4,5,5,5)
y<-c(1,2,4,1,2,3,1,3,1,1,2,3)
data<-c(3,2,1,6,4,1,2,3,1,1,1,3)

z<-cbind(x,y)

	m=max(z);
	if (is.array(z)) {
		n=length(dim(z)) 
		d=dim(z)
	} else {
		n=1
		d=length(data)
	}
	a=vector(length=m^n,mode="numeric")
	wts=1
	window=c(1,1)
	outDim=c(5,4)

s<-.C("kz_ma", as.double(a), as.double(outDim), as.double(data), as.integer(length(data)), as.double(z), as.integer(length(d)), as.double(d), as.double(window), NAOK=TRUE, DUP=FALSE, PACKAGE="kza")

/* Input :                                                          */
/*   outData  - empty double                                        */
/*   outDim  - output dimensions									*/
/*   inData - vector of data values									*/
/*   nIn - size of input vector										*/
/*   index - array of index locations								*/
/*   nDim - number of dimensions used by index						*/
/*	 inDim - array with size of each dimension of the index			*/
/*	 window - array with window filter dimensions					*/
/* data<-(1,2,3,4,5,6,7,8,9,10,11,12) */

/* 3 dimensions */
x<-c(1,1,1,2,2,2,3,3,4,5,5,5)
y<-c(1,2,4,1,2,3,1,3,1,1,2,3)
z<-c(2,3,1,2,2,3,1,2,1,1,2,3)
data<-c(3,2,1,6,4,1,2,3,1,1,1,3)

nz<-cbind(y,x,z)

	m=max(nz);
	if (is.array(nz)) {
		n=length(dim(nz)) 
		d=dim(nz)
	} else {
		n=1
		d=length(data)
	}
	a=vector(length=4*5*3,mode="numeric")
	wts=1
	window=c(0,0,0)
	outDim=c(4,5,3)
	nDim=d[2] 

s<-.C("kz_ma", as.double(a), as.double(outDim), as.double(data), as.integer(length(data)), as.double(nz), as.integer(nDim), as.double(d), as.double(window), NAOK=TRUE, DUP=FALSE, PACKAGE="kza")

b<-array(a, outDim)
dd<-melt(b)

dd<-subset(dd,!is.na(dd$value))
names(dd)<-c("y","x","z","data")

nz<-cbind(dd$y, dd$x, dd$z)
data<-dd$data


##################

x<-c(3,2,1,6,4,1,2,3,1,1,1,3)


x = as.vector(x) 
  n = length(x)
k = n
  k2 = k%/%2
  y=double(n)
 

.C("runmean_exact", x, y , as.integer(n), as.integer(k), NAOK=TRUE, DUP=FALSE, PACKAGE="kza") 