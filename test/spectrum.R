periodogram <- function(s) {
	fourier<-fft(s)
	
	magnitude<-Mod(fourier)
	 
	# extract the phase which is atan(Im(fourier)/Re(fourier))
	phase<-Arg(fourier)
	
	# select only first half of vectors
	magnitude_firsthalf <- magnitude[1:(length(magnitude)/2)]
	phase_firsthalf<-phase[1:(length(magnitude)/2)]
	 
	# generate x-axis
	x.axis <- 1:length(magnitude_firsthalf)/length(magnitude)
	
	p<-cbind(x.axis, magnitude_firsthalf)
	return(p)
}

#plot(p,type='l')