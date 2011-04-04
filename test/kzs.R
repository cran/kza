t <- seq(from = -round(400*pi), to = round(400*pi), by = .25) 

# Construct the signal over time
ts <- 0.5*sin(sqrt((2*pi*abs(t))/200))
signal <- ifelse(t < 0, -ts, ts)

# Bury the signal in noise [randomly, from N(0, 1)]
et <- rnorm(length(t), mean = 0, sd = 1)
yt <- et + signal

# Data frame of (t, yt) 
pts <- data.frame(cbind(t, yt))


# Apply 3 iterations of kzs
a<-system.time(kzs::kzs(y = pts[,2], x = pts[,1], smooth = 80, scale = .2, k = 3, edges = TRUE, plot = TRUE))
lines(signal ~ t, col = "red")
title(main = "kzs(smooth = 80, scale = .2, k = 3, edges = TRUE)")
legend("topright", c("True signal","kzs estimate"), cex = 0.8,
col = c("red", "black"), lty = 1:1, lwd = 2, bty = "n")

z<-kzft(yt,m=2,f=0,k=3,dim=1)
plot(Re(z),type='l')
lines(signal, col="red")

b<-kzs(yt)

################
## example 2
################
obs <- seq(1:length(t))
t20 <- sample(obs, size = length(obs)/5)
pts20 <- pts[-t20,]        

# Plot of (t,yt) with 20 percent of the data removed
plot(pts20$yt ~ pts20$t, main = "Signal buried in noise\n20 percent of 
(t, yt) deleted", xlab = "t", ylab = "yt", type = "p")

# Apply 3 iterations of kzs
kzs::kzs(y = pts20[,2], x = pts20[,1], smooth = 80, scale = .2, k = 3, edges = TRUE, plot = TRUE)
lines(signal ~ t, col = "red")
title(main = "kzs(smooth = 80, scale = .2, k = 3, edges = TRUE)")
legend("topright", c("True signal","kzs estimate"), cex = 0.8, 
col = c("red", "black"), lty = 1:1, lwd = 2, bty = "n")  



function (y, x, smooth, scale, k = 1, edges = TRUE, plot = TRUE) 
{
    x <- sort(x)
    dx <- diff(x)
    if (scale > min(dx[dx > 0])) 
        stop("'scale' must be less than or equal to the minimum difference of consecutive X values")
    if (smooth >= (max(x) - min(x))) 
        stop("'smooth' must be much less than the difference of the max and min X values")
    h <- smooth/2
    xrange <- range(x)
    for (i in 1:k) {
        xi <- as.vector(x)
        maxx <- max(xi)
        minx <- min(xi)
        yvals <- y
        xk <- seq(minx - h, maxx + h, scale)
        yk <- numeric(length(xk))
        for (j in 1:length(xk)) {
            w <- abs(xi - xk[j])
            Ik <- which(w <= h)
            size <- length(Ik)
            Yik <- yvals[Ik]
            yk[j] <- (1/size) * sum(Yik)
        }
        df <- data.frame(cbind(xk, yk))
        data <- na.omit(df)
        x <- data$xk
        y <- data$yk
    }
    if (edges == FALSE) {
        edgs <- data[data$xk >= min(xrange) & data$xk <= max(xrange), 
            ]
        data <- as.data.frame(edgs)
    }
    if (plot == TRUE) {
        plot(data$yk ~ data$xk, xlab = "xk", ylab = "yk", xlim = c(min(xrange), 
            max(xrange)), type = "l", pch = 19, col = "black")
    }
    return(data)
}




t2<-1:3000
f<-10/1000
x<-cos(2*pi*f*t2) + rnorm(length(t2),0,2)

z3<-kzft(x,m=2,f=0, k=3, dim=1)
z3<-kzft(x,m=2,f=0, k=3, dim=1)


