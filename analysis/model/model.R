### The main purpose of this file is to run and compare different models (that are usually nested).
### Since other datasets will lead the analyst down different paths.
### This is difficult to understand if you're not familar with running mamximum likelihood and Bayesian
### multilevel models in R and OpennBUGS.  Much of these decisions are informed by the second half of Gelman & Hill (2007).

rm(list=ls(all=TRUE))
library(arm) #Load Gelman & Hill's multilevel model package.
library(R2OpenBUGS) #Communicate with OpenBUGS from R
library(MCMCpack) #For 'rwish' function.

pathData  <- "./data/raw/qpcr.tsv"
pathModel <- "./analysis/model/model.txt"

ds <- read.table(pathData, header=TRUE, sep="\t")
#summary(ds)
colnames(ds) <- c("Treatment", "PlotID", "Extraction",  "Month", "Year", "Bacteria", "LogCopyCount", "SoilMoisture", "SoilMoisturePlot")
ds <- subset(ds, PlotID != 1) #Drop Plot #1
ds$Year <- ds$Year + 2003 #Reexpress so the first year is 2004
ds$Month <- 12 - ds$Month * 4 #Reexpress so month '2' is April and month '1' is August
useMedian <- FALSE
ds$Bacteria <-relevel(ds$Bacteria, "Total") #Make 'Total' the first level for the later contrasts.
ds$SoilMoisturePlot[is.na(ds$SoilMoisturePlot)] <- 10 #Assume all mising moisture measurements were because it was too low to measure; assign a low value.

ds <- cbind(ds, Time=ds$Year + ds$Month/12, LogCopyCountResidual=rep(NA,times=nrow(ds)))
#ds <- cbind(ds, CenteredMoisture=ds$SoilMoisturePlot - mean(ds$SoilMoisturePlot, na.rm=T))
#ds <- cbind(ds, LogCenteredMoisture=log(ds$CenteredMoisture))
ds <- cbind(ds, LogMoisture=log(ds$SoilMoisturePlot))
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
#ds[ds$Treatment==1 & ds$Year==2005 & ds$Month==8 & ds$PlotID %in% c(3,4) & ds$Extraction==2 & ds$Bacteria=='Actinobacteria' & ds$LogCopyCountResidual < -.9, 'LogCopyCount'] <- NA

ds$PlotID <- factor(ds$PlotID)
ds$Extraction <- factor(ds$Extraction)
ds$Treatment <- factor(ds$Treatment)
ds$Month <- factor(ds$Month)
ds$Year <- factor(ds$Year)
ds$Time <- factor(ds$Time)
ds <- subset(ds, !is.na(LogCopyCount))

#(m1 <- lmer(LogCopyCount ~ Bacteria + CenteredMoisture +Treatment +(0 + CenteredMoisture | PlotID) + (0 + Treatment | Time ), data=ds))
#(m2 <- lmer(LogCopyCount ~ Bacteria + CenteredMoisture +Treatment +(1 + CenteredMoisture | PlotID) + (0 + Treatment | Time ), data=ds))
#(m17b <- lmer(LogCopyCount ~ Bacteria + CenteredMoisture +Treatment +(1 + CenteredMoisture | PlotID) + (1 + Treatment | Time ), data=ds))
#(m21 <- lmer(LogCopyCount ~ Bacteria + CenteredLogMoisture +Treatment +(1 + CenteredLogMoisture | PlotID) + (1 + Treatment | Time ), data=ds))
#ranef(m21)
#se.ranef(m21)
#(m22 <- lmer(LogCopyCount ~ Bacteria + CenteredLogMoisture +Treatment +(1  | PlotID) + (1 + Treatment | Time ), data=ds))
#ranef(m22)
#se.ranef(m22)
#(m23 <- lmer(LogCopyCount ~ Bacteria + Treatment +(1  | PlotID) + (1 + Treatment | Time ), data=ds))
#ranef(m23)
#(m24 <- lmer(LogCopyCount ~ Bacteria + CenteredLogMoisture +Treatment +(1 | PlotID) + (1 + Treatment + CenteredLogMoisture  | Time ), data=ds))
#ranef(m24)
#(m25 <- lmer(LogCopyCount ~ Bacteria + TimeCenteredLogMoisture +Treatment +(1 | PlotID) + (1 + Treatment  | Time ), data=ds))
#ranef(m25)
#(m26 <- lmer(LogCopyCount ~ Bacteria + TimeCenteredLogMoisture +Treatment +(1 | PlotID) + (1 + Treatment + TimeCenteredLogMoisture  | Time ), data=ds))
#ranef(m26)
#(m27 <- lmer(LogCopyCount ~ Bacteria + TimeCenteredLogMoisture + Treatment +(1 | PlotID) + (1 + Treatment | Time)+ (0 + TimeCenteredLogMoisture | Time), data=ds))
#ranef(m27)
#(m <- lmer(LogCopyCount ~ 1 + ( 1 | PlotID ), data=ds))

#m28 <- lmer(LogCopyCount ~ CenteredLogMoisture +Treatment + (1 | Bacteria) + (1 | PlotID) + (1 + Treatment | Time ), data=ds)
#m28 #Current favorite
#ranef(m28)
#
#(m29 <- lmer(LogCopyCount ~ CenteredLogMoisture +Treatment + (1+CenteredLogMoisture | Bacteria) + (1 +CenteredLogMoisture | PlotID) + (1 +CenteredLogMoisture+ Treatment | Time ), data=ds))
#ranef(m29)

bugsFit <- NA ; rm(bugsFit); (startTime <- Sys.time())
selectionVector <- !is.na(ds$CenteredLogMoisture) #& ds$Bacteria=="Total"
y <- ds$LogCopyCount[selectionVector]
n <- length(y)
plotID <- as.numeric(ds$PlotID[selectionVector])
plotCount <- length(unique(plotID))
time <- as.numeric(ds$Time[selectionVector])
timeCount <- length(unique(time))
bac <- as.numeric(ds$Bacteria[selectionVector])
bacCount <- length(unique(bac))
moisture <- ds$CenteredLogMoisture[selectionVector]
tx <- as.numeric(ds$Treatment[selectionVector]) - 1
W <- diag(2)
data <- list("y", "n", "plotID", "plotCount", "time", "timeCount", "bac", "bacCount", "tx", "W", "moisture")
inits <- function( ) { list(
   sigma.y=runif(1), mu0=rnorm(1, mean=7), mu1=rnorm(1), mu2=rnorm(1)
  , bacMean=rnorm(bacCount)#, rho1=runif(1, -.97, .97)
  #, B=array(rnorm(2*plotCount), c(plotCount,2))
  , D=array(rnorm(2*timeCount), c(timeCount,2))
  , Tau.D=rwish(3, diag(2))
  , sigma.b0=runif(1)#, sigma.b1=runif(1)
  , xi.d0=runif(1, min=.11, max=5), xi.d2=runif(1, min=.11, max=5)
)}
parameters <- c( "mu0","mu0.adj", "sigma.y", "mu1", "mu2", "bacOffset"
  ,"b0.adj", "sigma.b0"#, "b1", "sigma.b1", "rho1" #, "b0"
  , "d0.adj", "d2", "rho2" #, "sigma.d0", "sigma.d1", "d0"
  , "xi.d0", "xi.d2", "sigma.d0.adj", "sigma.d2.adj"
  , "txNotOffset"#, "droughtDifTx","droughtDifControl"
  )
(bugsFit <- bugs(
  data                   = data,
  inits                  = inits,
  parameters.to.save     = parameters,
  model.file             = pathModel,
  n.chains               = 6,
  n.iter                 = 80, #8000*11
  n.burnin               = 2, #2000
  n.thin                 = 10,
  working.directory      = getwd(),
  debug                  = TRUE
))
#(bugsFit <- bugs(data, inits, parameters, pathModel, n.chains=6, n.iter=1000, bugs.directory=pathBugs, debug=T))

# x6 took 6.81218 hours
# x7 took 8.486831 hours with competing background tasks.
#n.thin=10 & fixed moisture:
# x6 took 7.3 hours (lots of background tasks)
# x8 took 7.15

Sys.time() - startTime

#The follow features are important conclusions of the paper.  They were all true for almost all the models,
#   which suggests the patterns/features aren't artifacts of a particular model, and are likely to be replicated.
#mean(bugsFit$sims.list$d2[,1]<0)   #Feature A
#mean(bugsFit$sims.list$d2[,3]>0)   #Feature B
#mean(bugsFit$sims.list$droughtDifTx > 0 ) #Feature C
#mean(bugsFit$sims.list$droughtDifControl < 0 ) #Feature D
#mean(bugsFit$sims.list$featureCProb)
max(bugsFit$summary[,'Rhat'])
