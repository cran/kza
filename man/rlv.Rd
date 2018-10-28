\name{rlv}
\alias{rlv}
\title{Rolling Variance Function}
\description{
  Rolling variance function, (rlv) will detect data variations of 1-dimensional
  vector, 2-dimensional matrix and/or 3-dimensional array, which in turn display
  their local discontinuities (boundaries).
}
\details{
	The rlv is calculated based on a dynamic rolling window. Within a given
	window, directly estimating its variance allows to capture the local
	variations of the signal underlying analysis. Since this window is moving,
	when its size is optimal, the significant local variance atlas will be
	captured. It is neither necessary to make any assumption, nor to model
	anything, so it is a nonparametric statistic. However, direct application of
	rlv to data embedded in heavy background noises could fail to identify the
	significant local variations. In situations when the estimated data are
	embedded in background noise, necessary smoothness may be expected.
}
\usage{
rlv(inpt, krnl)
}
\arguments{
  \item{inpt}{input data.}
  \item{krnl}{length of the rolling window for detecting.}
}
\references{
  Igor G Zurbenko, Mingzeng Sun. Estimation of spatial boundaries with rolling
  variance and 2D KZA algorithm. Biometrics & Biostatistics International
  Journal. 7(4):263.270, 2018
}
\examples{
## A. Signal/objects
# to create a cylinder
circleFun <- function(center = c(0,0),diameter = 1, npoints = 200){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints) 
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
rxo1 <- circleFun(center = c(38, 38), diameter = 38, npoints = 200)
rxi1 <- circleFun(center = c(38, 38), diameter = 18, npoints = 200)
rxo2 <- circleFun(center = c(58, 58), diameter = 38, npoints = 200)
rxi2 <- circleFun(center = c(58, 58), diameter = 18, npoints = 200)
rxo3 <- circleFun(center = c(28, 78), diameter = 18, npoints = 200)
rxi3 <- circleFun(center = c(28, 78), diameter = 8, npoints = 200)

set.seed(3)
s2 <- matrix(rep(0,100*100),nrow=100) 
for(i in 1:nrow(rxo1)) {
s2[rxo1[i,2],38:rxo1[i,1]]=0.5
}
for(i in 1:nrow(rxi1)) {
s2[rxi1[i,2],38:rxi1[i,1]]=1
}
# second cylinder
for(i in 1:nrow(rxo2)) {
s2[rxo2[i,2],58:rxo2[i,1]]=0.5
}
for(i in 1:nrow(rxi2)) {
s2[rxi2[i,2],58:rxi2[i,1]]=1
}
#third separated cylinder
for(i in 1:nrow(rxo3)) {
s2[rxo3[i,2],28:rxo3[i,1]]=0.3
}
for(i in 1:nrow(rxi3)) {
s2[rxi3[i,2],28:rxi3[i,1]]=0.6
}
signal=s2
## Graphing Library plotly
\dontrun{
plot_ly(z=signal, type="surface")
}

## B. Rolling variance atlas

rolv=rlv(s2,3)
\dontrun{
plot_ly(z=rolv, type="surface")
}

}
\keyword{rlv}
