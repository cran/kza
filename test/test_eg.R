#####
# wedge example (3d)
# This can take quite a bit of time to run
####
w<-array(data=0, dim=c(100,100,100))
for (i in 35:65) {
    w[i:65,35:65,i] <- 1
}
w<-w+rnorm(n = 100*100*100, sd = 2)
m=20
system.time(a<-kz(w,m,3))
system.time(b<-kza(w,y=a,m=m,k=3,impute_tails=TRUE))

######################
# movie of filtered object
#####################
x11()
for(i in 1:50) {
    image(matrix(b$kza[,,i],50,50),col=gray(seq(0,min(max(b$kza[,,i],na.rm=TRUE),1),1/255)))
}



########################
# example using index set
########################
obs <- seq(1:length(t))
t20 <- sample(obs, size = length(obs)/5)
pts20 <- pts[-t20,]        

b<-kzs(yt,m=400)
par(mfrow=c(2,1))
plot(yt,type='l')
lines(signal,col="red")

plot(b,type='l', ylim=c(-0.5,1))
lines(signal,col="red")
par(mfrow=c(1,1))
