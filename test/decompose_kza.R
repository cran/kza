sdecomp_kza<-function(y, k=3, period=c(53,4))
{
    if (!is.ts(y)) stop("Time series object required.")
    
    l <- length(y)
    f <- frequency(y)
    m=period[1]
    
    trend<-ts(kza(y,q=m,k=k),start=start(y), freq=f)
    i=which.max(diff(trend))
    cp<-i
	b=max(diff(trend))
	a<-y[1:i]
	a2<-y[(i+1):length(y)]-b

    seasonal_freq=period[2]
    
    #s1<-kz(append(a,a2),q=3,k=1)
	l<-length(a)
    season<-ts(kz(a-trend[1:l],q=seasonal_freq,k=k),start=start(y),frequency=f)

    periods <- l%/%f
    index <- c(cumsum(rep(f, periods - 1)))
    figure <- numeric(f)
    for (i in 1L:f) figure[i] <- mean(season[index + i])   
    
    one_season = ts(figure,start=start(y), frequency=f)
    seasonal_mean <- ts(rep(figure, periods), start = start(y), frequency = f)
	seasonal_mean <- ts(kz(seasonal_mean,q=2,k=3), start = start(y), frequency = f)
 	s<-start(seasonal_mean)[2]
 	j<-f-s+2
	one_season<-seasonal_mean[j:(j+f)]
	
    trend2<-ts(kz(y,q=m,k=k),start=start(y), freq=f)
    season<-ts(kz(y-trend2,q=seasonal_freq,k=k),start=start(y),frequency=f)
    fitted<-ts(season+trend2, start=start(y), frequency=f)

    structure(list(fitted = fitted,
            x = y,
            trend = trend,
            seasonal_trend = season,
            seasonal_mean = seasonal_mean,
            one_season = one_season,
            frequency = seasonal_freq,
            change_point = cp,
            residuals = y-fitted
            ),
        class = "kza")
}

predict.kza <-
    function (object, n.ahead = 1, level = 0.95, ...)
{
    f <- frequency(object$x)
    trend<-object$trend
    fitted<-object$fitted
    
    b<-object$change_point
    a<-ts(trend[(b+1):length(trend)],start=time(trend)[(b+1)],frequency=frequency(trend))
    
	a2<-lm(a~time(a))
	p<-predict(a2,n.ahead=n.ahead)

	d<-p[1]-a[length(a)]
	p=p[1:n.ahead]-d
	
    et<-end(trend)[2]
    es<-length(object$one_season)-1
    
    s1<-object$one_season[(et+1):(es)]
    ee<-es-length(s1)
    s2<-object$one_season[1:ee]

	ns<-append(s1,s2)
	
	figure<-ns
	periods <- n.ahead%/%f
    for (i in 1L:periods) figure <- append(figure,ns)
    
    s<-figure[1:n.ahead]
    fit<-p+s
	
    return (fit)    
}

plot.kza <-
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

plot.decompose.kza <- function(x, ...)
{
    plot(cbind(
               observed=x$x,
               trend    = x$trend,
               seasonal = x$seasonal_trend,
               fitted   = x$fitted
               ),
         main = paste("KZ Decomposition of time series"),
         ...)
}

