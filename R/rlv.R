
##Authors: Mingzeng Sun	 Igor G Zurbenko
##Department of Epidemiology and Biostatistics, State University of New York at Albany, New York, USA
##Email: msun@albany.edu; igorg.zurbenko@gmail.com

##August 18, 2018 reviewed
##it is required that krnl is odd number 
rlv=function(inpt, krnl){
inpt_chk=ncol(rbind(dim(inpt)))
if (length(inpt_chk)==0) {
# 1D vector
#note: make the working space with 2*krnl bigger than the data space at all directions
l=length(inpt)+2*krnl
imgn0= rep(0,l)
imgn1=imgn0
imgn1[(krnl+1):(l-krnl)]=inpt

rolv=imgn0
for (i in (krnl+1):(length(imgn0)-krnl)) {

ssp0 <- imgn1[(i-(0.5*(krnl-1))):(i+(0.5*(krnl-1)))]
ssp0a <- imgn1[(i-1):i] ## left, including self
ssp0b <- imgn1[i:(i+1)] ## right, including self
if (krnl>1) {vrn0=sd(ssp0)*sd(ssp0)  
}
if (krnl==1) {vrn0a=sd(ssp0a)*sd(ssp0a)  
vrn0b=sd(ssp0b)*sd(ssp0b)  
vrn0=max(vrn0a, vrn0b)  
}
rolv[i]=vrn0
}
rolVariance=rolv[(krnl+1):(length(rolv)-krnl)]
} else if (length(inpt_chk)==1 & inpt_chk==2) {
# 2D matrix
#note: make the working space with 2*krnl bigger than the data space at both directions
l=nrow(inpt)+2*krnl
w=ncol(inpt)+2*krnl
imgn0= matrix(rep(0,l*w),nrow=l)
imgn1=imgn0
imgn1[(krnl+1):(l-krnl), (krnl+1):(w-krnl)]=inpt

rolv=imgn0
# starting loop
for (i in (krnl+1):(nrow(imgn0)-krnl)) {
for (j in (krnl+1):(ncol(imgn0)-krnl)) {
ssp0 <- imgn1[(i-(0.5*(krnl-1))):(i+(0.5*(krnl-1))), (j-(0.5*(krnl-1))):(j+(0.5*(krnl-1)))]
ssp0a <- imgn1[(i-1):i, j:(j+1)] ## right-above corner, including self
ssp0b <- imgn1[(i-1):i, (j-1):j] ## left-above corner, including self
ssp0c <- imgn1[i:(i+1), (j-1):j] ## left-bottom corner, including self
ssp0d <- imgn1[i:(i+1), j:(j+1)] ## right-bottom corner, including self

if (krnl>1) {vrn0=sd(ssp0)*sd(ssp0)  
}
if (krnl==1) {vrn0a=sd(ssp0a)*sd(ssp0a)  
vrn0b=sd(ssp0b)*sd(ssp0b)  
vrn0c=sd(ssp0c)*sd(ssp0c)  
vrn0d=sd(ssp0d)*sd(ssp0d)
vrn0=max(vrn0a, vrn0b, vrn0c, vrn0d)  
}
rolv[i,j]=vrn0
}
}
rolVariance=rolv[(krnl+1):(nrow(rolv)-krnl), (krnl+1):(ncol(rolv)-krnl)]

} else if (length(inpt_chk)==1 & inpt_chk==3) {
# 3D array
##it is required that krnl is odd number 
#note: make the working space with 2*krnl bigger than the data space at both directions
l=nrow(inpt)+2*krnl
w=ncol(inpt)+2*krnl
h=2*krnl+length(inpt)/(nrow(inpt)*ncol(inpt))
imgn0= array(rep(0,l*w*h),dim=c(l, w, h))
imgn1=imgn0
imgn1[(krnl+1):(l-krnl), (krnl+1):(w-krnl), (krnl+1):(h-krnl)]=inpt

rolv=imgn0
# starting loop 1
for (i in (krnl+1):(l-krnl)) {
for (j in (krnl+1):(w-krnl)) {
for (k in (krnl+1):(h-krnl)) {
ssp0 <- imgn1[(i-(0.5*(krnl-1))):(i+(0.5*(krnl-1))), (j-(0.5*(krnl-1))):(j+(0.5*(krnl-1))), (k-(0.5*(krnl-1))):(k+(0.5*(krnl-1)))]
ssp0a <- imgn1[(i-1):i, (j-1):j, (k-1):k] ## one of the 8 corner, including self
ssp0b <- imgn1[(i-1):i, (j-1):j, k:(k+1)] ## one of the 8 corner, including self
ssp0c <- imgn1[(i-1):i, j:(j+1), (k-1):k] ## one of the 8 corner, including self
ssp0d <- imgn1[(i-1):i, j:(j+1), k:(k+1)] ## one of the 8 corner, including self
ssp0e <- imgn1[i:(i+1), (j-1):j, (k-1):k] ## one of the 8 corner, including self
ssp0f <- imgn1[i:(i+1), (j-1):j, k:(k+1)] ## one of the 8 corner, including self
ssp0g <- imgn1[i:(i+1), j:(j+1), (k-1):k] ## one of the 8 corner, including self
ssp0h <- imgn1[i:(i+1), j:(j+1), k:(k+1)] ## one of the 8 corner, including self

if (krnl>1) {vrn0=sd(ssp0)*sd(ssp0)  
}
if (krnl==1) {vrn0a=sd(ssp0a)*sd(ssp0a)  
vrn0b=sd(ssp0b)*sd(ssp0b)  
vrn0c=sd(ssp0c)*sd(ssp0c)  
vrn0d=sd(ssp0d)*sd(ssp0d)
vrn0e=sd(ssp0e)*sd(ssp0e)
vrn0f=sd(ssp0f)*sd(ssp0f)
vrn0g=sd(ssp0g)*sd(ssp0g)
vrn0h=sd(ssp0h)*sd(ssp0h)
vrn0=max(vrn0a, vrn0b, vrn0c, vrn0d, vrn0e, vrn0f, vrn0g, vrn0h)  
}
rolv[i,j,k]=vrn0
}
}
}
rolVariance=rolv[(krnl+1):(l-krnl), (krnl+1):(w-krnl), (krnl+1):(h-krnl)]
}
return(rolVariance)
}


