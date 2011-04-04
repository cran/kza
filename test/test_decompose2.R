decompose.kz<-function(y, k=3, period=c(365))
{
    if (!is.ts(y)) stop("Time series object required.")
    
    l <- length(y)
    f <- frequency(y)
    m=period[1]

    trend<-kz(y,q=m,k=k)
    base_trend<-ts(trend[(2*f):(length(trend)-2*f-1)],start=start(y),frequency=f)

    s<-start(base_trend)[1]+start(base_trend)[2]/f
    e<-end(base_trend)[1]+end(base_trend)[2]/f
    t<-seq(s,by=1/f,length=length(base_trend))
    model<-lm(base_trend~t)
    t<-append(t,seq(e+1/f,by=1/f,length=2*f))
    p<-predict(model,data.frame(x=I(t)))
    p2<-p[(length(base_trend)+1):length(p)]
    slide<-p2[1]-base_trend[length(base_trend)]
    p2<-p2-slide

    t<-seq(s,by=1/f,length=length(base_trend))
    model<-lm(rev(base_trend)~t)
    t<-append(t,seq(e+1/f,by=1/f,length=2*f))
    p<-predict(model,data.frame(x=I(t)))
    p1<-p[(length(base_trend)+1):length(p)]
    slide<-p1[1]-base_trend[1]
    p1<-rev(p1-slide)

    trend<-ts(c(p1,base_trend,p2),start=start(y), frequency=f)
    
    a<-kzp(y-trend,m)
    s<-which.max(a)
    seasonal_freq=(s-1)/m
    
    season<-ts(2*Re(kzft(y-trend,m=m,k=3,dim=1,trim=)),start=start(y),frequency=f)

    periods <- l%/%f
    index <- c(cumsum(rep(f, periods - 3)))
    figure <- numeric(f)
    for (i in 1L:f) figure[i] <- mean(season[index + i])   
    
    one_season = ts(figure,start=start(y), frequency=f)
    seasonal_mean <- ts(rep(figure, periods), start = start(y), frequency = f)
    
    season<-append(one_season,season)
    season<-append(season,one_season)

    fitted<-ts(season+trend, start=start(y), frequency=f)

    structure(list(fitted = fitted,
            x = y,
            trend = trend,
            seasonal_trend = season,
            seasonal_mean = seasonal_mean,
            one_season = one_season,
            frequency=seasonal_freq
            ),
        class = "kz")
}

predict.kz <-
    function (object, n.ahead = 1, level = 0.95, ...)
{
    f <- frequency(object$x)
    trend<-object$trend
    
    s<-start(trend)[1]+start(trend)[2]/f
    e<-end(trend)[1]+end(trend)[2]/f
    t<-seq(s,by=1/f,length=length(trend))
    model<-lm(trend~t)
    t<-append(t,seq(e+1/f,by=1/f,length=n.ahead))
    p<-predict(model,data.frame(x=I(t)))
    p<-p[(length(trend)+1):length(p)]
    
    start=1
    s<-object$one_season[start:n.ahead]
    fit<-p+s

    return (fit)    
}
