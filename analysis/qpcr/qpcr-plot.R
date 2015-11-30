#library(bootstrap)
#rm(list=ls(all=TRUE))


xRange <- c(2004.64, 2006.35)
yRange <- c(-.5,.5) #Used for aggregates of plots -with bands (hinges)
moistureResidualScale <- 0.1
boxColor <- gray(.8)

TreatmentColorsSolid <- function( index ) { return( c('blue', 'tomato')[index] ) }
TreatmentColorsAlpha <- function( index ) {
  transparency <- .15
  return( c(rgb(0, 0, 1, alpha=transparency), rgb(1, .3882353, .2784314, alpha=transparency))[index] )
}

#Max Rhat 1.018272
#Wishart, Fixed Moisture
controlUpper <-c(-.05347, -.004079, .3365, -.05325)
controlIntercepts <- c(-.1149, -.05472, .2854, -.1158)
controlLower <-c(-.1764, -.1052, .2338, -.1784)
txUpper <- c(.1972, -.1158, -.176, -.1601)
txIntercepts <- c(.1256, -.1855, -.2398, -.2473)
txLower <- c(.0538, -.2553, -.304, -.3345)

##Wishart, Random Moisture
#controlUpper <-c(-.03434, -.04281, .3506, -.05475)
#controlIntercepts <- c(-.09519, -.09523, .3031, -.1127)
#controlLower <-c(-.1571, -.1432, .2557, -.1739)
#txUpper <- c(.2224, -.1333, -.161, -.1947)
#txIntercepts <- c(.1544, -.2049, -.2163, -.2819)
#txLower <- c( .08377, -.2776, -.2786, -.3768)
#
#NonWishart
#controlUpper <-c(-.03122, -.04554, .3514, -.05131)
#controlIntercepts <- c(-.09499, -.09612, .3031, -.122)
#controlLower <-c(-.1576, -.1479, .2528, -.1748)
#txUpper <- c(.2214, -.1332, -.1562, -.1921)
#txIntercepts <- c(.1546, -.2046, -.2166, -.2849)
#txLower <- c( .08965, -.2766, -.276, -.3708)

#txUpper <- c(.1467, 2.319, -.2227, 2.789)
#txIntercepts <- c(.04813, -.35, -.311, -.4785)
#txLower <- c(-.04303, -3.048, -.4014, -2.91)
##  fixedTxSe <- 0 #se.fixef(model)['Treatment2']   #Should this be included?
#  controlSEs <- sqrt(fixedTxSe^2 + se.ranef(model)$Time[, 1]^2)
#  txSEs <- sqrt(fixedTxSe^2 + se.ranef(model)$Time[, 2]^2)
#
op <- par(mfrow=c(1,1), mar=c(3.5,2,2.5,2.5))
plot(NA, xlim=xRange, ylim=yRange, bty="n", xaxt="n", yaxt="n", xaxs="i", mgp=c(1.1,0.1,0),
  xlab="", ylab="log 16s rRNA Gene Copies Residuals",
  main="MLM Expected Residuals", sub="95% MLM CIs around Tx mean for each time")

yTicks <- c(-.5, 0, .5)
axis(2, at=c(-.5, 0, .5), mgp=c(1.1,0.2,0), tcl=.25)
axis(2, at=c(-.5,-.25, 0,.25, .5), mgp=c(1.1,0.2,0), tcl=.15, labels=rep("",5))
xValues <- c(2004.667, 2005.333, 2005.667, 2006.333)
axis(1, at=xValues, labels=c("Aug\n2004","Apr\n2005","Aug\n2005","Apr\n2006"), tck=0)
tickLocations <- seq(from=2005, to=2006, by=1)
axis(1, at=tickLocations, rep(NA,length=length(tickLocations)))
box(col=gray(.8))   

lines(x=xValues, y=controlIntercepts, col=TreatmentColorsSolid(1), lwd=2) 
lines(x=xValues, y=txIntercepts, col=TreatmentColorsSolid(2), lwd=2)
points(x=xValues, y=controlIntercepts, col=TreatmentColorsSolid(1), cex=.9, pch=16)
points(x=xValues, y=txIntercepts, col=TreatmentColorsSolid(2), cex=.9, pch=16)

bandControl <- c(controlUpper,  rev(controlLower))
bandTx <- c(txUpper,  rev(txLower))

polygon(x=c(xValues, rev(xValues)), y=bandControl, col=TreatmentColorsAlpha(1), border=NA)
polygon(x=c(xValues, rev(xValues)), y=bandTx, col=TreatmentColorsAlpha(2), border=NA)



par(op)
