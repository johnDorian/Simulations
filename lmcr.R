#!/usr/bin/Rscript
################################################################################
################################################################################
### Title: LMCR_FINAL.R ########################################################
### Author: Jason Lessels j.lessels@usyd.edu.au ################################
### Created: Wed Feb  3 15:42:52 2010 ##########################################
### Modified: Wed Feb  3 15:42:52 2010 #########################################
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
det.fit <- function(object){
tpn <- object$model$TP$psill[1]
fln <- object$model$FLOW$psill[1]
tpfln <- object$model$TP.FLOW$psill[1]
d <- matrix(c(tpn,tpfln,tpfln,fln),c(2,2))
print(det(d))
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
data$TP <- data$TP+0.01
data$FLOW <- data$FLOW+0.01
## Determine the correct lambda values for the transformations of both TP and Flow
TP.lambda <- boxcox.fit(data$TP)
FLOW.lambda <- boxcox.fit(data$FLOW)
data <- data
data$TP <- bct(data$TP,TP.lambda$lambda)
data$FLOW <- bct(data$FLOW,FLOW.lambda$lambda)

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
cross<-v[1:bins,3]
auto.1<-v[(bins+1):(2*bins),3]
auto.2<-v[(2*bins+1):(3*bins),3]
dist<-v[(bins+1):(2*bins),2]
pairs<-v[(bins+1):(2*bins),1]
semvar=data.frame(dist,auto.1,cross,auto.2,pairs)

#lmcr<-function(semvar,nolags,nvar,wgt,icvp,cparf,modtyp,covar,maxdist,guessa,lock)
#need to know wgt,cparf,covar,guessa,lock
lmcr(semvar,20,2,?,1,?,4,?,200,?,?)

#Plot the cross covariogram and get some initial values for it.
plot(v,g)
#The cross covariance models using initial values, these will change.
g = gstat(g,id=c("TP","FLOW"),model=vgm(0.5,"Sph",48,0.2),fill.all=TRUE)

### STEP 3. Find the range of each covariogram using fit.lmc
##Use fit.lmc() function to find the ranges for each covariogram.
fit.range=fit.lmc(v,g,fit.lmc=TRUE,fit.ranges=TRUE)
fit.range
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
