######
#######
## f=12 monthly
## f=52 weekly
## f=365 or 366 daily
## f=365*24 hourly

###############
## daily signal 
###############
decompose.kz<-function(y, k=3, freq=c(365,30))
{
	if (!is.ts(y)) stop("Time series object required.")
	
	l <- length(y)
	f <- frequency(y)
	m=freq[1]

	trend<-kz(y,q=m,k=k)
	base_trend<-ts(trend[(2*f):(length(trend)-2*f-1)],start=start(y),frequency=f)

	model<-HoltWinters(base_trend)
	rmodel<-HoltWinters(ts(rev(base_trend),start=start(y),frequency=f))
	
	a2<-predict(model,n.ahead=(2*f))
	a1<-rev(predict(rmodel,n.ahead=(2*f)))
	trend<-ts(c(a1,base_trend,a2),start=start(y), frequency=f)
	
	# monthly
	seasonal_freq=1/freq[2]
	if (f==365 || f==366) s<-freq[2]
	season<-ts(2*Re(kzft(y-trend,m=s,k=k,f=seasonal_freq,dim=1)),start=start(y),frequency=f)
	#season<-ts(2*Re(kz(y-trend,q=s,k=k)),start=start(y),frequency=f)

	periods <- l%/%f
    index <- c(cumsum(rep(f, periods - 3)))
    figure <- numeric(f)
    for (i in 1L:f) figure[i] <- mean(season[index + i])   
    
	one_season = ts(figure,start=start(y), frequency=f)
	seasonal_mean <- ts(rep(figure, periods), start = start(y), frequency = f)

	fitted<-ts(season+trend, start=start(y), frequency=f)

	structure(list(fitted = fitted,
			x = y,
			trend = trend,
			seasonal_trend = season,
			seasonal_mean = seasonal_mean,
			one_season = one_season,
			random = ts(y-fitted,start=start(y),frequency=frequency(y))
			),
		class = "kz")
}

predict.kz <-
    function (object, n.ahead = 1, level = 0.95, ...)
{
 	f <- frequency(object$x)
	trend<-object$trend

	model<-HoltWinters(trend)
	p<-predict(model,n.ahead=n.ahead)
	
	start=1
	s<-object$one_season[start:n.ahead]
	fit<-p+s

	return (fit)	
}

plot.kz <-
    function (x, predicted.values = NA, intervals = TRUE, separator = TRUE,
              col = 1, col.predicted = 2, col.intervals = 4, col.separator = 1,
              lty = 1, lty.predicted = 1, lty.intervals = 1, lty.separator = 3,
              ylab = "Observed / Fitted", main = "Kolmogorov-Zurbenko filtering",
              ylim = NULL, ...)
{
    if (is.null(ylim))
      ylim <- range(na.omit(c(fitted(x), x$x, predicted.values)))

    preds <- length(predicted.values) > 1 || !is.na(predicted.values)

    ## plot fitted/predicted values
    plot(ts(c(fitted(x), if(preds) predicted.values[,1]),
            start = start(fitted(x)),
            frequency = frequency(fitted(x))),
         col = col.predicted,
         ylim = ylim,
         ylab = ylab, main = main,
         lty = lty.predicted,
         ...
         )

    ## plot observed values
    lines(x$x, col = col, lty = lty)
}

plot.decompose.kz <- function(x, ...)
{
    plot(cbind(
               observed=x$x,
               trend    = x$trend,
               seasonal = x$seasonal_mean,
               random   = x$random
               ),
         main = paste("KZ Decomposition of time series"),
         ...)
}

