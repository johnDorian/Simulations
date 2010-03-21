#!/usr/bin/Rscript
################################################################################
################################################################################
### Title: LMCR_FINAL.R ########################################################
### Author: Jason Lessels j.lessels@usyd.edu.au ################################
### Created: Wed Feb  3 15:42:52 2010 ##########################################
### Modified: Wed Feb  21/03/2010 ##############################################
### Aim: To perform LMCR #######################################################
### Description: This script performs LMCR on temporal data ####################
################################################################################
### NOTES: This still isn't great but it's the best there is. A further note is
### in the original script (lmcr_cont.R) likfit was used to determine the vario-
### gram variables, but this was found not to make any difference in the end re-
### sults. Therefore variofit is still for the cross-variogram only. ###########
################################################################################
### TODO: Run through all the possiblities #####################################
################################################################################

###GET SOME USEFUL FUNCTIONS OUT OF THE WAY.
## Let's universalise plotting on all Oprating Systems.
win.graph = function(){if(.Platform$OS=='unix') x11() else win.graph=win.graph}
## Now so we can easly check the determinates of fit.lmc using det.fit
covStructure <- function(object){
auto1 <- object$model[1][[1]][[2]][1]
auto2 <- object$model[3][[1]][[2]][1]
cross <- object$model[2][[1]][[2]][1]
d <- matrix(c(auto1,cross,cross,auto2),c(2,2))
return(d)
}

###Now we can actually start with the script.

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
## Add arbituary value to the data, so we can transform the data.
data$TP <- log(data$TP+0.01)
data$FLOW <- log(data$FLOW+0.01)

##Create gstat objects of the data
spdf <- SpatialPointsDataFrame(data[,1:2],data)

### STEP 1. Create A gstat object with TP and flow.
g = gstat(NULL, "TP", spdf$TP ~ 1, spdf,,maxdist=200)
g = gstat(g, "FLOW", spdf$FLOW ~ 1, spdf,maxdist=200)

### STEP 2. Create intial values for the models.
#Create variogarms for plotting.
v = variogram(g,cutoff=200,width=200/20)
####	Need to re-organise for the lmcr function.
bins<-as.numeric(summary(v$id))[1]
cross<-split(v,v$id)[[1]]$gamma
auto.1<-split(v,v$id)[[2]]$gamma
auto.2<-split(v,v$id)[[3]]$gamma
dist<-split(v,v$id)[[1]]$dist
pairs<-split(v,v$id)[[1]]$np
semvar=data.frame(dist,auto.1,cross,auto.2,pairs)

#lmcr<-function(semvar,nolags,nvar,wgt,icvp,cparf,modtyp,covar,maxdist,guessa,lock)
#need to know wgt,cparf,covar,guessa,lock


#Plot the cross covariogram and get some initial values for it.
plot(v,g)
#The cross covariance models using initial values, these will change.
g = gstat(g,id=c("TP","FLOW"),model=vgm(0.5,"Exp",48,0.2),fill.all=TRUE)

### STEP 3. Find the range of each covariogram using fit.lmc
##Use fit.lmc() function to find the ranges for each covariogram.
fit=fit.lmc(v,g,fit.lmc=FALSE,fit.range=TRUE)
covar<-covStructure(fit)
eigen(covar)
lmcr(semvar,20,2,1,1,cparf,4,covar,200,48,0)



##Avergae the ranges out and use the average range as the max range.
range<-mean(fit.range$model$TP[[3]][2],fit.range$model$TP.FLOW[[3]][2],fit.range$model$FLOW[[3]][2])
g.ave.range = gstat(g,id=c("TP","FLOW"),model=vgm(0.5,"Sph",range,0.2),fill.all=TRUE)
##Create another model using 2 days (48) as the range).
g.2days.range=gstat(g,id=c("TP","FLOW"),model=vgm(0.5,"Sph",48,0.2),fill.all=TRUE)

### STEP 4. Fit the two final models using fit.lmc()
fit.ave=fit.lmc(v,g.ave.range,fit.lmc=TRUE,fit.ranges=FALSE)
fit.2days=fit.lmc(v,g.2days.range,fit.lmc=TRUE,fit.ranges=FALSE)
##Check the determinates of the model make sure they are positive
det.fit(fit.ave)
det.fit(fit.2days)
##Plotting the results.
pdf("lmc_with_ave_range.pdf")
plot(v,fit.ave,main="Average range used for model")
dev.off()
win.graph()
pdf("lmc_with_2day_range.pdf")
plot(v,fit.2days,main="2 day range used for model")
##Compare the values for the four best models.
fit.ave
fit.2days
