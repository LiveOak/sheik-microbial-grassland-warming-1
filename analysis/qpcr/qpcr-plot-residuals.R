### The main purpose of this file is to inspect the models that are estimated invidually for
### each bacterial group.  If the groups resemble each other, then it *may* make sense to unify
### them in the same multilevel model.  When looking for similarities, pay attention to
### (a) relative locations of the mean & median points within a tx condition,
### (b) relative locations of the mean & median points between tx conditions, and
### (c) spread/variance between tx conditions at a given time point.

### The MLM panel should have much smaller spread than the first seven panels, because
### it is pulling from all seven groups (among other reasons).

rm(list=ls(all=TRUE))
library(bootstrap) #Load Efron & Tibshirani's bootstrap package.
library(arm) #Load Gelman & Hill's multilevel model package.

path <- "./data/raw/qpcr.tsv"
bootstrapSize <- 1000 #Use 1,000 for debugging; use 10,000 for publication

ds <- read.table(path, header=TRUE, sep="\t")
#summary(ds)
colnames(ds) <- c("Treatment", "PlotID", "Extraction",  "Month", "Year", "Bacteria", "LogCopyCount", "SoilMoisture", "SoilMoisturePlot")
ds <- subset(ds, PlotID!=1) #Drop Plot #1
ds$Year <- ds$Year + 2003 #Reexpress so the first year is 2004
ds$Month <- 12 - ds$Month * 4 #Reexpress so that month '2' is April and month '1' is August
useMedian <- FALSE
ds$Bacteria <-relevel(ds$Bacteria, "Total") #Make 'Total' the first level for the later contrasts.
ds$SoilMoisturePlot[is.na(ds$SoilMoisturePlot)] <- 10 #Assume all mising moisture measurements were because it was too low to measure; assign a low value.

ds <- cbind(ds, Time=ds$Year + ds$Month/12, LogCopyCountResidual=rep(NA,times=nrow(ds)))
ds <- cbind(ds, LogMoisture=log(ds$SoilMoisture))
ds <- cbind(ds, CenteredLogMoisture=ds$LogMoisture - mean(ds$LogMoisture, na.rm=T))
timePointCount <- length(unique(ds$Time))

for( bacteriaName in unique(ds$Bacteria) ) {
  for( plotIndex in unique(ds$PlotID) ) {
    dsPlot <- subset(ds, ds$PlotID==plotIndex & ds$Bacteria==bacteriaName & !is.na(LogCopyCount))
    if( useMedian ) centerLogCopy <- median(dsPlot$LogCopyCount, na.rm=TRUE)
    else centerLogCopy <- mean(dsPlot$LogCopyCount, na.rm=TRUE)
    ds$LogCopyCountResidual[ds$PlotID==plotIndex & ds$Bacteria==bacteriaName] <- ds$LogCopyCount[ds$PlotID==plotIndex & ds$Bacteria==bacteriaName] - centerLogCopy
  }
}
#Filter out two abonormally low values
ds[ds$Treatment==1 & ds$Year==2005 & ds$Month==8 & ds$PlotID %in% c(3,4) & ds$Extraction==2 & ds$Bacteria=='Actinobacteria' & ds$LogCopyCountResidual < -.9, c('LogCopyCount', 'LogCopyCountResidual')]<- rep(NA, 2)

xRange <- c(2004.63, 2006.36) #range(ds$Time)
yRange <- c(-.5,.5) #Used for aggregates of plots -with bands (hinges)
moistureResidualScale <- 0.1
boxColor <- gray(.8)

PlotColors <- function( index ) { return( rainbow(length(unique(ds$PlotID)))[index] ) }
TreatmentColorsSolid <- function( index ) { return( c('blue', 'tomato')[index] ) }
TreatmentColorsAlpha <- function( index ) {
  transparency <- .15
  return( c(rgb(0, 0, 1, alpha=transparency), rgb(1, .3882353, .2784314, alpha=transparency))[index] )
}

plotColors <- rainbow(length(unique(ds$PlotID))) #plotColors[2] <- "brown"
treatmentDot <- c(16, 3) #A plus sign for tx and a small dot for control.
treatmentLinePatterns <- c(1, 1) #treatmentLinePatterns <- c(2, 1)
treatmentErorBarsOffset <- c(0, 0) #treatmentErorBarsOffset <- c(-.02, .02)

PlotLabelsAndTicks <- function( panelLabel, displayXLabels, displayYAxisValues ) {
  yTicks <- c(-.4, 0, .4)
  if( displayXLabels ) {
    xValuesLabels <- c(2004.7, 2005.333, 2005.667, 2006.25)
    axis(1, at=xValuesLabels, labels=c("Aug\n2004","Apr\n2005","Aug\n2005","Apr\n2006"), tck=0)
    tickLocations <- seq(from=2005, to=2006, by=1)
    axis(1, at=tickLocations, rep(NA,length=length(tickLocations)))
  }
  if( displayYAxisValues ) {
    axis(2, at=yTicks, mgp=c(1.1,0.1,0), tcl=.25, col.axis=gray(.6))
  }

  #Plot the left and right vertical axes
  axis(2, at=c(-.4,-.2, 0,.2, .4), mgp=c(1.1,0.2,0), tcl=.15, labels=rep("",5), col=gray(.8)) #Short unlabeled gray ticks
  axis(2, at=c(yTicks), mgp=c(1.1,0.2,0), tcl=.3, labels=rep("",3)) #longer labeled black ticks
  axis(4, at=c(-.4,-.2, 0,.2, .5), mgp=c(1.1,0.2,0), tcl=.15, labels=rep("",5), col=gray(.8))
  axis(4, at=c(yTicks), mgp=c(1.1,0.2,0), tcl=.3, labels=rep("",3))
  box(col=boxColor)

  text(x=2004.66, y=.44, panelLabel, pos=4, xpd=NA )#Annotate each panel with the bacterial group.
}

#This function is used for the first seven panels (ie, for each bacterial group).
PlotBacteria <- function( bacteria, panelLabel, displayYAxisValues=T, displayYLabel=FALSE, displayXLabels=FALSE ) {
  plot(NA, xlim=xRange, ylim=yRange, bty="n", xaxt="n", yaxt="n", xaxs="i",
    mgp=c(1.1,0.1,0), xlab="", ylab="", main="")
  PlotLabelsAndTicks(panelLabel, displayXLabels, displayYAxisValues)

  for( treatmentIndex in unique(ds$Treatment) ) {
    longitudinalXValues <- rep(NA, length=length(unique(ds$Time)))
    longitudinalYValues <- rep(NA, length=length(unique(ds$Time)))
    longitudinalYValuesSoil <- rep(NA, length=length(unique(ds$Time)))
    upperBand <- rep(NA, length=length(unique(ds$Time)))
    lowerBand <- rep(NA, length=length(unique(ds$Time)))
    longitudinalIndex <- 1

    for( timeIndex in unique(ds$Time) ){ #?Sort time
      dsSlice <- subset(ds, Treatment==treatmentIndex & Time==timeIndex  & !is.na(LogCopyCountResidual) & Bacteria==bacteria)
      if(useMedian) yValues <- median(dsSlice$LogCopyCountResidual)
      else yValues <- mean(dsSlice$LogCopyCountResidual)
      xValues <- rep(timeIndex, length=length(yValues)) + treatmentErorBarsOffset[treatmentIndex] #+ rnorm(n=length(yValues), sd=.015)
      points(x=xValues, y=yValues, pch=treatmentDot[treatmentIndex], col=TreatmentColorsSolid(treatmentIndex), lwd=2)

      longitudinalXValues[longitudinalIndex] <- xValues
      longitudinalYValues[longitudinalIndex] <- yValues

      #Inspect the spread as defined by the (a)hinges, (b) parametric CI, and (c) bootstrap CI
      #Hinges
      upperBand[longitudinalIndex] <- quantile(dsSlice$LogCopyCountResidual, probs=.75)
      lowerBand[longitudinalIndex] <- quantile(dsSlice$LogCopyCountResidual, probs=.25)

      #Standard parametric CIs
      standardError <- sqrt(var(dsSlice$LogCopyCountResidual) / length(dsSlice$LogCopyCountResidual))
      upperBand[longitudinalIndex] <- yValues + 2 * standardError
      lowerBand[longitudinalIndex] <- yValues - 2 * standardError

      #Bootstrap CIs
      #results <- bcanon(dsSlice$LogCopyCountResidual, nboot=bootstrapSize, mean, alpha=c(.025, .975))
      #upperBand[longitudinalIndex] <- results$confpoints[2, 2]
      #lowerBand[longitudinalIndex] <- results$confpoints[1, 2]
      longitudinalIndex <- longitudinalIndex + 1
    }#End Time loop

    lines(x=longitudinalXValues, y=longitudinalYValues,  lty=treatmentLinePatterns[treatmentIndex], col=TreatmentColorsSolid(treatmentIndex), lwd=2)
    vertices <- c(upperBand, rev(lowerBand))
    polygon(x=c(unique(ds$Time), rev(unique(ds$Time))), y=vertices, col=TreatmentColorsAlpha(treatmentIndex), border=NA)
  }#End Treatment
} #End of function

#This function is used for the last panel (ie, the mean and spread of the predicted model).
PlotMlm <- function( bacteriaType, displayYLabel=F ) {
  plot(NA, xlim=xRange, ylim=yRange, bty="n", xaxt="n", yaxt="n", xaxs="i", mgp=c(1.1,0.1,0),
    xlab="", ylab="", main="", sub="95% MLM CIs around Tx mean for each time")
  PlotLabelsAndTicks("H) MLM Model", displayXLabels=T, displayYAxisValues=F)

  #These coefficients came from the Bayesian model that had a Wishart prior and a fixed effect for Moisture.
  #These are the vertical coordinates for the intercepts and the CI bands (ie, upper & lower).
  controlUpper <-c(-.039, -.004261, .3382, -.05393)
  controlIntercepts <- c(-.1151, -.0545, .2856, -.1159)
  controlLower <-c(-.1783, -.1048, .2328, -.1781)
  txUpper <- c(.1973, -.1176, -.1746, -.1592)
  txIntercepts <- c(.1253, -.186, -.2399, -.2478)
  txLower <- c(.05314, -.2555, -.304, -.3344)

  xValues <- c(2004.667, 2005.333, 2005.667, 2006.333)
  #Draw the points and connected them with lines.
  lines(x=xValues, y=controlIntercepts, col=TreatmentColorsSolid(1), lwd=2)
  lines(x=xValues, y=txIntercepts, col=TreatmentColorsSolid(2), lwd=2)
  points(x=xValues, y=controlIntercepts, pch=treatmentDot[1], col=TreatmentColorsSolid(1), lwd=2)
  points(x=xValues, y=txIntercepts, pch=treatmentDot[2], col=TreatmentColorsSolid(2), lwd=2)

  #Draw the CI bands
  polygon(x=c(xValues, rev(xValues)), y=c(controlUpper, rev(controlLower)), col=TreatmentColorsAlpha(1), border=NA)
  polygon(x=c(xValues, rev(xValues)), y=c(txUpper, rev(txLower)), col=TreatmentColorsAlpha(2), border=NA)
}

oldPar <- par(mfrow=c(2,4), mar=c(0,1,2.2,0.2))  #The top four panels don't have a bottom margin.
PlotBacteria(bacteria='Total', panelLabel='A) Total Bacteria')
PlotBacteria(bacteria='Actinobacteria', panelLabel='B) Actinobacteria', displayYAxisValues=F)
PlotBacteria(bacteria='Acidobacteria', panelLabel='C) Acidobacteria')
PlotBacteria(bacteria='Alphaproteobacteria', panelLabel='D) Alphaproteobacteria', displayYAxisValues=F)

par(mar=c(2.2,1,0,0.2))#The bottom four panels don't have a top margin.
PlotBacteria(bacteria='Crenarcheota', panelLabel='E) Crenarcheota', displayXLabel=T, displayYLabel=T)
PlotBacteria(bacteria='Planctomycetes', panelLabel='F) Planctomycetes', displayXLabel=T, displayYAxisValues=F)
PlotBacteria(bacteria='Verrucomicrobia', panelLabel='G) Verrucomicrobia', displayXLabel=T)
PlotMlm( bacteriaType="All Types")
par(oldPar) #Reset the graphical parameters.

#If you're not running a Bayesian model, convert these variables to factors for the MLM routine in the 'lmer' package.
#But make sure they're not factors when you graph, because they will no longer be numeric variables.
#ds$PlotID <- factor(ds$PlotID)
#ds$Extraction <- factor(ds$Extraction)
#ds$Treatment <- factor(ds$Treatment)
#ds$Month <- factor(ds$Month)
#ds$Year <- factor(ds$Year)
#ds$Time <- factor(ds$Time)
#ds <- subset(ds, !is.na(LogCopyCount))



