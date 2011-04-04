
kz_decompose<-function(y, k=3, period=NULL)
{
    if (!is.ts(y)) stop("Time series object required.")

    l <- length(y)
    f <- frequency(y)
    
	# if period is not given determine periods for daily or weekly frequency
	if (is.null(period)) {
		if (f==365 || f==366) period=c(365,8) # filter out weekly signal
		if (f==52 || f==53) period=c(f,4)
	}

    m=period[1]
    trend<-ts(kz(y,m=m,k=k),start=start(y), freq=f)

    seasonal_freq=period[2]
    season<-ts(kz(y-trend,m=seasonal_freq,k=k),start=start(y),frequency=f)

    periods <- l%/%f
    index <- c(cumsum(rep(f, periods - 1)))
    figure <- numeric(f)
    for (i in 1L:f) figure[i] <- mean(season[index + i])   
    
    one_season = ts(figure,start=start(y), frequency=f)
    seasonal_mean <- ts(rep(figure, periods), start = start(y), frequency = f)
	seasonal_mean <- ts(kz(seasonal_mean,m=2,k=3), start = start(y), frequency = f)
 	s<-start(seasonal_mean)[2]
 	j<-f-s+2
	one_season<-seasonal_mean[j:(j+f)]

    fitted<-ts(season+trend, start=start(y), frequency=f)

    structure(list(
    		observed=y, seasonal=season, seasonal_mean=seasonal_mean, trend=trend, fitted=fitted, residuals=y-fitted,
            one_season = one_season,
            seasonal_frequency=seasonal_freq,
            call=match.call()
            ),
        class = "kzd")
}

#alias
kzd<-kz_decompose

predict.kzd <-
    function (object, n.ahead = 1, prediction.interval = FALSE, level = 0.95, ...)
{
    f <- frequency.kzd(object)
    trend<-trend.kzd(object)
    fitted<-fitted.kzd(object)
    
    hw<-HoltWinters(trend, gamma=FALSE)

	p<-predict(hw,n.ahead=n.ahead, level=level, prediction.interval=prediction.interval)
#	d<-p[1]-trend[length(trend)]
#   p=p[1]-d
    
    t1<-end(trend)[2]
    t2<-n.ahead+t1-1
    
    s<-object$one_season
    iterations = trunc(t2/length(object$one_season))
    while (iterations>0) {
    	s<-append(s,s)
    	iterations=iterations-1
	}    	
    	
    fit<-p[,1]+s[t1:t2]
	
    return (fit)    
}

plot.kzd <- function(x, ...)
{
    plot(cbind(
               observed = cbind(observed=x$observed, 
               		seasonal=x$seasonal, 
               		seasonal_mean=x$seasonal_mean, 
               		trend=x$trend, 
               		fitted=x$fitted, 
               		residuals=x$residuals)
               ),
         main = paste("KZ Decomposition of time series"), 
         ...)
}

summary.kzd <- function(object, digits = getOption("digits"), ...)
{
    cat(" Call:\n ")
    dput(object$call, control=NULL)
    cat("\n Time.series components:\n")
    print(summary(cbind(object$observed, object$seasonal, object$seasonal_mean, object$trend, object$fitted, object$residuals), digits = digits, ...))

    invisible(object)
}

residuals.kzd <- function(object, ...) 
{
	return (object$residuals)
}

fitted.kzd <- function(object, ...) 
{
	return (object$fitted)
}

trend.kzd <- function(object, ...) 
{
	return (object$trend)
}

frequency.kzd <- function(x, ...) 
{
	return (x$seasonal_frequency)
}

mean_kz <- function(x, frequency, scale, offset=0, periods=NULL)
{
    if (is.null(periods)) {periods <- length(x)%/%frequency}
    index <- c(0, cumsum(rep(frequency, periods - 1))-1)
    figure <- numeric(frequency)
    for (i in 1L:frequency) figure[i] <- mean(x[index + i + offset])   
    one_season = ts(figure, start=start(x), frequency=scale)
}
