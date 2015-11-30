### The main purpose of this file is to plot soil moisture and bacterial abundance in the same space.
### No coefficients are estimated here (except for the CIs around moisture).  That happens in the other files.

rm(list=ls(all=TRUE))
library(bootstrap) #Load Efron & Tibshirani's bootstrap package.

dateRange <- as.numeric(as.Date(c('2004/2/1', '2006/5/1'))) #Define the minimum and maximum
path <- "./data/raw/soil-moisture.csv"
ds <- read.table(path, header=TRUE, sep="\t")
colnames(ds) <- c('Date', 'PlotID', 'Rep', 'Treatment', 'Moisture')
ds <- subset(ds, PlotID %in% 2:5)
useColor <- TRUE #Toggle the color on/off for presentations/articles.
bootstrapSize <- 1000 #Use 1,000 for debugging; use 10,000 for publication

ds$Date <- as.Date(ds$Date, format='%m/%d/%Y')
ds <- subset(ds, dateRange[1] <= Date & Date <= dateRange[2])
ds$Date <- as.numeric(ds$Date)
ds <- ds[order(ds$Date), ] #Make sure they're sorted by date, for the graph
#head(ds, 10) #Uncomment to inspect the first ten lines of the dataset.
#summary(ds) #Uncomment to inspect the variables of the dataset.

uniqueDays <- unique(ds$Date, na.rm=T)
uniqueDays <- uniqueDays[!is.na(uniqueDays)]
boxColor <- gray(.8)
labelColor <- gray(.6)

bacteriaDates <- as.numeric(as.Date(c('2004/08/15', '2005/04/15', '2005/08/15', '2006/04/15')))

#Residuals: bacteriaControl <- c(-0.1963724,3.859266e-02, 0.2633013, 0.0873106)
#Residuals: bacteriaTreatment <- c( 5.269522e-02, -0.03250291, -0.1991696, -0.01385490)
bacteriaControl <- c(-.1151, -.0545, .2856, -.1159) #Bayes w/ fixed moisture (Coefficients are estimated with WinBUGS)
bacteriaTreatment <- c(.1253, -.186, -.2399, -.2478) #Bayes w/ fixed moisture (Coefficients are estimated with WinBUGS)

if( useColor ) {
  TreatmentColorsSolid <- function( index ) { return( c('blue', 'tomato')[index] ) }
  TreatmentColorsAlpha <- function( index ) {
    transparency <- .15
    return( c(rgb(0, 0, 1, alpha=transparency), rgb(1, .3882353, .2784314, alpha=transparency))[index] )
  }
  TreatmentColorsAlpha2 <- function( index ) {
    transparency <- .3
    return( c(rgb(0, 0, 1, alpha=transparency), rgb(1, .3882353, .2784314, alpha=transparency))[index] )
  }
}
if( !useColor ) {
  TreatmentColorsSolid <- function( index ) { return( c('black', 'gray')[index] ) }
  TreatmentColorsAlpha <- function( index ) {
    transparency <- .15
    return( c(rgb(.5, .5, .5, alpha=transparency), rgb(.5, .5, .5, alpha=transparency))[index] )
  }
  TreatmentColorsAlpha2 <- function( index ) {
    transparency <- .3
    return( c(rgb(0, 0, 0, alpha=transparency), rgb(.745098, .745098, .745098, alpha=transparency))[index] )
  }
}
treatmentDot <- c(16, 3)
treatmentErorBarsOffset <- c(0, 0)#To nudge them horizontally use something like c(-.02, .02)
yRange <- range(ds$Moisture, na.rm=T)

secondaryDateLocations <- as.Date(c("2004-04-15", "2004-12-15", "2005-12-15"))
secondaryDateLabels <- c("Apr\n2004", "Dec\n2004", "Dec\n2006")

#plot(ds$Date, ds$Moisture) #Very crude longitudinal graph of moisture (make sure nothing screwy happened with the fancy graph.)
oldPar <- par(mfrow=c(1,1), mar=c(1.9,1.5,0,1.7), family="sans")
plot(NA, xlim=dateRange, ylim=yRange, bty="n", xaxt="n",yaxt="n", xaxs="i", yaxs="i",#yaxt="n",
  mgp=c(1,0.1,0), tcl=.25, cex.axis=.8,cex.lab=.8, main="", xlab="", ylab="")

yTicks <- c(10, 20, 30, 40)
axis(2, at=yTicks, mgp=c(1,0,0), tcl=.25, col=boxColor, cex.axis=.8) #Plot vertical axis
mtext(side=2, "Percent Moisture", line=.7)

dateLabels <- as.Date(c('2004/08/15', '2005/04/15', '2005/08/15', '2006/04/15'))
axis(1, at=as.numeric(dateLabels), labels=c("Aug\n2004","Apr\n2005","Aug\n2005","Apr\n2006"), tck=0, cex.axis=.8, padj=-.17)
axis(1, at=secondaryDateLocations, labels=secondaryDateLabels, tck=0, col.axis=boxColor, cex.axis=.8, padj=-.17)
dateTicks <- unique(ds$Date)
axis(1, at=as.numeric(dateTicks), labels=rep(NA,length=length(dateTicks)) , tcl=.25, col=boxColor)
box(col=boxColor)

for( treatmentIndex in unique(ds$Treatment) ) {#Cycle through both treatment conditions
  longitudinalXValues <- rep(NA, length=length(unique(ds$Time)))
  longitudinalYValues <- rep(NA, length=length(unique(ds$Time)))
  longitudinalYValuesSoil <- rep(NA, length=length(unique(ds$Time)))
  upperBand <- rep(NA, length=length(unique(ds$Time)))
  lowerBand <- rep(NA, length=length(unique(ds$Time)))
  longitudinalIndex <- 1

  for( timeIndex in uniqueDays ){ #Cycle through the days (make sure time is sorted above)
    dsSlice <- subset(ds, Treatment==treatmentIndex & Date==timeIndex  & !is.na(Moisture) ) #Select a subset of the data.
    yValues <- dsSlice$Moisture
    xValues <- rep(timeIndex, length=length(yValues))#+ treatmentErorBarsOffset[treatmentIndex]
    #Uncomment line below to see the raw data
    #points(x=xValues, y=yValues, pch=treatmentDot[treatmentIndex], col=TreatmentColorsAlpha(treatmentIndex))

    longitudinalXValues[longitudinalIndex] <- xValues[1]
    longitudinalYValues[longitudinalIndex] <- mean(dsSlice$Moisture, na.rm=T)

    #Inspect the spread as defined by the (a)hinges, (b) parametric CI, and (c) bootstrap CI
    #Hinges
    upperBand[longitudinalIndex] <- quantile(dsSlice$Moisture, probs=.75)
    lowerBand[longitudinalIndex] <- quantile(dsSlice$Moisture, probs=.25)

    #Standard parametric CIs
    #standardError <- sqrt(var(dsSlice$Moisture) / length(dsSlice$Moisture))
    #upperBand[longitudinalIndex] <-  longitudinalYValues[longitudinalIndex] + 2 * standardError
    #lowerBand[longitudinalIndex] <-  longitudinalYValues[longitudinalIndex] - 2 * standardError

    #Bootstrap CIs
    #results <- bcanon(dsSlice$Moisture, nboot=bootstrapSize, mean, alpha=c(.025, .975))
    #upperBand[longitudinalIndex] <- results$confpoints[2, 2]
    #lowerBand[longitudinalIndex] <- results$confpoints[1, 2]

    longitudinalIndex <- longitudinalIndex + 1
  }#End Time loop

  #Plot the mean moisture points and connect them with lines.
  points(x=longitudinalXValues, y=longitudinalYValues, pch=treatmentDot[treatmentIndex], col=TreatmentColorsSolid(treatmentIndex), lwd=2)
  lines(x=longitudinalXValues, y=longitudinalYValues, col=TreatmentColorsSolid(treatmentIndex), lwd=2)

  #Construct the upper and lower boundaries for the translucent CI band.
  vertices <- c(upperBand, rev(lowerBand))
  polygon(x=c(unique(ds$Date), rev(unique(ds$Date))), y=vertices, col=TreatmentColorsAlpha(treatmentIndex), border=NA)
}#End Treatment loop

#This function scales the bacteria abundance values so they're visible on a plot showing the moisture data.
#   In other words, choose 'scaleCoefficient' so the dashed lines (and the right axis) fits.
RescaleBacteriaValues <- function( abundance ) {
  scaleCoefficient <- 60
  return( scaleCoefficient * abundance + median(ds$Moisture, na.rm=T) )
}

#Plot the points for Bayesian estimated model values and connect them with lines.
points(x=bacteriaDates, y=RescaleBacteriaValues(bacteriaControl), col=TreatmentColorsAlpha2(1), pch=treatmentDot[1], lwd=2, lty="42")
points(x=bacteriaDates, y=RescaleBacteriaValues(bacteriaTreatment), col=TreatmentColorsAlpha2(2), pch=treatmentDot[2], lwd=2, lty="F4")
lines(x=bacteriaDates, y=RescaleBacteriaValues(bacteriaControl), col=TreatmentColorsAlpha2(1), lwd=2, lty="44")
lines(x=bacteriaDates, y=RescaleBacteriaValues(bacteriaTreatment), col=TreatmentColorsAlpha2(2), lwd=2, lty="FF")

#Draw and label the right side axis
bacteriaTicks <- c(-.5, -.25, 0, .25, .5) #
axis(4, at=RescaleBacteriaValues(bacteriaTicks), labels=bacteriaTicks, mgp=c(.9,-.3, 0), tcl=.25, col=boxColor, col.axis=boxColor, cex.axis=.8)
mtext("Abundance of Total Bacterial (MLM Residuals)", side=4, line=.5, cex=.8, col=labelColor)

#Annotate the graph with the text "Drought Begins"
if( useColor ) colorDroughtBegins <- "brown" else colorDroughtBegins <- "black"
mtext("Drought\nBegins", side=c(1,3), at=as.numeric(as.Date("2005-04-28")), line=-.1, col=colorDroughtBegins, las=2, adj=c(0, 1), cex=.8)

par(oldPar) #Reset the graphical parameters.
