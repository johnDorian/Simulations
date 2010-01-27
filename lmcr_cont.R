#!/usr/bin/Rscript
################################################################################
################################################################################
### Title: LMCR_CONT.R #########################################################
### Author: Jason Lessels j.lessels@usyd.edu.au ################################
### Created: Wed Jan 27 14:15:45 2010 ##########################################
### Modified: Wed Jan 27 14:15:45 2010 #########################################
### Aim: To use the geoR and the gstat packages to fit a variogram #############
### Description: This script will use geoR-lifit fn. for guidence. #############
################################################################################
### NOTES: This script will re-use the lmcr.R script but add more to it. #######
################################################################################
### TODO: Have to make the script from various scripts in the past. ############
################################################################################

## Let's universalise plotting on all Oprating Systems.
win.graph = function(){if(.Platform$OS=='unix') x11() else win.graph()=x11()}
## Set the working directory.
setwd("~/Documents/code/Simulations/")
## Load the datasets.
data<-read.csv("TP_flow.csv", header=TRUE)
## Load the nessacry libraries. A couple for variograms, and some for the boxcox transformations.
library(gstat);library(geoR);library(MASS);library(TeachingDemos);library(moments)
## Clean up the dataset - create a back.up of the data
back.up <- data
data <- data[,-c(1,2)]
names(data) <- c("X","Y","TP","FLOW")
data$TP <- data$TP+0.01
data$FLOW <- data$FLOW+0.01
## Determine the correct lambda values for the transformations of both TP and Flow
TP.lambda <- boxcox.fit(data$TP)
FLOW.lambda <- boxcox.fit(data$FLOW)
## See if the lambda values are any good.
win.graph()
par(mfrow=c(2,2))
hist(log(data$TP),xlab="skewness = 0.96",main="Log Transformation of TP")
hist(log(data$FLOW),xlab="skewness = -0.45",main="Log Transformation of Flow")
hist(bct(data$TP,TP.lambda$lambda),xlab="skewness = 0.067, lambda = -0.42",main="Boxcox Trans. of TP")
hist(bct((data$FLOW),FLOW.lambda$lambda),xlab="skewness = 0.016, lambda = -0.088",main="Boxcox Trans. of TP")
## Get some values for the transformations
skewness(data$TP)
skewness(data$FLOW)
skewness(log(data$TP))
skewness(log(data$FLOW))
skewness(bct(data$TP,TP.lambda$lambda))
skewness(bct(data$FLOW,FLOW.lambda$lambda))

## Now we need to use the geoR package to perform likfit for both TP and Flow using both transformations.
## Create two datasets for each transformation method.
data.log <- data
data.log$TP <- log(data$TP)
data.log$FLOW <- log(data$FLOW)
data.bc <- data
data.bc$TP <- bct(data$TP,TP.lambda$lambda)
data.bc$FLOW <- bct(data$FLOW,FLOW.lambda$lambda)
## Create 4 geodata objects one four each variable and transformation.
qual.tp.log<-as.geodata(data.1,coords.col=1:2,data.col=3)
qual.flow.log<-as.geodata(data.1,coords.col=1:2,data.col=4)
qual.tp.bc<-as.geodata(data.2,coords.col=1:2,data.col=3)
qual.flow.bc<-as.geodata(data.2,coords.col=1:2,data.col=4)
## Perform likfit on each different variable
log.tp.likfit <- likfit(qual.tp.log,ini=c(0.5,0.5),lik.method="REML")
log.flow.likfit <- likfit(qual.flow.log,ini=c(0.5,0.5),lik.method="REML")
bc.tp.likfit <- likfit(qual.tp.bc,ini=c(0.5,0.5),lik.method="REML")
bc.flow.likfit <- likfit(qual.flow.bc,ini=c(0.5,0.5),lik.method="REML")
## Have a look at the summaries of all the likfit models.
summary(log.tp.likfit)
summary(log.tp.likfit)
summary(bc.flow.likfit)
summary(bc.flow.likfit)

## Rip the psill from the likfits
log.tp.psill <- summary(log.tp.likfit)$spatial.component[[2]][1]
log.flow.psill <- summary(log.flow.likfit)$spatial.component[[2]][1]
bc.tp.psill <- summary(bc.tp.likfit)$spatial.component[[2]][1]
bc.flow.psill <- summary(bc.flow.likfit)$spatial.component[[2]][1]

## Rip the psill from the likfits
log.tp.nug <- summary(log.tp.likfit)$nugget.component[[2]]
log.flow.nug <- summary(log.flow.likfit)$nugget.component[[2]]
bc.tp.nug <- summary(bc.tp.likfit)$nugget.component[[2]]
bc.flow.nug <- summary(bc.flow.likfit)$nugget.component[[2]]


## Rip the range values from the likfits
log.tp.range <- summary(log.tp.likfit)$practicalRange
log.flow.range <- summary(log.flow.likfit)$practicalRange
bc.tp.range <- summary(bc.tp.likfit)$practicalRange
bc.flow.range <- summary(bc.flow.likfit)$practicalRange


## Use the above results of the likfit to run the lmcr code.
log.spdf <- SpatialPointsDataFrame(data.log[,1:2],data.log)
bc.spdf <- SpatialPointsDataFrame(data.bc[,1:2],data.bc)
#Calculate the experimental auto-semivariogram
log.tp.vgm <- variogram(log.spdf$TP~1,log.spdf)
log.flow.vgm <- variogram(log.spdf$FLOW~1,log.spdf)
bc.tp.vgm <- variogram(bc.spdf$TP~1,bc.spdf)
bc.flow.vgm <- variogram(bc.spdf$FLOW~1,bc.spdf)

## 
win.graph()
par(mfrow=c(2,2))
plot(log.tp.vgm[,2],log.tp.vgm[,3])
plot(log.flow.vgm[,2],log.flow.vgm[,3])
plot(bc.tp.vgm[,2],bc.tp.vgm[,3])
plot(bc.flow.vgm[,2],bc.flow.vgm[,3])

##
log.tp.fit <- fit.variogram(log.tp.vgm,model=vgm(log.tp.psill,"Sph",log.tp.range,log.tp.nug))
log.flow.fit <- fit.variogram(log.flow.vgm,model=vgm(log.flow.psill,"Sph",log.flow.range,log.flow.nug))
bc.tp.fit <- fit.variogram(bc.tp.vgm,model=vgm(bc.tp.psill,"Sph",bc.flow.range,bc.tp.nug))
bc.flow.fit <- fit.variogram(bc.flow.vgm,model=vgm(bc.flow.psill,"Sph",bc.flow.range,bc.flow.nug))

win.graph()
par(mfrow=c(2,2))

plot(log.tp.vgm,log.tp.fit)
plot(log.flow.vgm,log.flow.fit)
plot(bc.tp.vgm,bc.tp.fit)
plot(bc.flow.vgm,bc.flow.fit)

###MODELLING LINEAR MODEL OF COERGIONALISATION
## Completed in two sections - one for log and the second for boxcox.
## The first one is for logged data.

#Creates a dataset for auto-semivariogram estimation
g.log = gstat(NULL, "TP", qual.sp.1$TP ~ 1, qual.sp.1)
#Somehow links b1 with b0 so it knows cross-semivariograms are needed
g.log = gstat(g.log, "FLOW", qual.sp.1$FLOW ~ 1, qual.sp.1)
#Calculates auto and cross-semivariograms
v = variogram(g.log,cutoff=200,width=200/20)
plot(v)
#Defines form of LMCR - range, variogram model
g.log = gstat(g.log, model = vgm(1, "Sph", 3000, 1), fill.all = TRUE)
#Fits LMCR model
g.fit = fit.lmc(v, g.log)

plot(v, g.fit)

#Creates a dataset for auto-semivariogram estimation
g.bc = gstat(NULL, "TP", qual.sp.2$TP ~ 1, qual.sp.2)
#Somehow links b1 with b0 so it knows cross-semivariograms are needed
g.bc = gstat(g.log, "FLOW", qual.sp.2$FLOW ~ 1, qual.sp.2)
#Calculates auto and cross-semivariograms
v = variogram(g.bc,cutoff=200,width=200/20)
plot(v)
#Defines form of LMCR - range, variogram model
g.bc = gstat(g.bc, model = vgm(1, "Sph", 3000, 1), fill.all = TRUE)
#Fits LMCR model
g.fit = fit.lmc(v, g.bc)

plot(v, g.fit)
