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
## Transform the data using two methods (natural log) and boxcox
data.log <- data
data.log$TP <- log(data$TP)
data.log$FLOW <- log(data$FLOW)
## Determine the correct lambda values for the transformations of both TP and Flow
TP.lambda <- boxcox.fit(data$TP)
FLOW.lambda <- boxcox.fit(data$FLOW)
data.bc <- data
data.bc$TP <- bct(data$TP,TP.lambda$lambda)
data.bc$FLOW <- bct(data$FLOW,FLOW.lambda$lambda)

##Create gstat objects of the data
log.spdf <- SpatialPointsDataFrame(data.log[,1:2],data.log)
bc.spdf <- SpatialPointsDataFrame(data.bc[,1:2],data.bc)
## Create two gstat objects for each transformation
g.log = gstat(NULL, "TP", log.spdf$TP ~ 1, log.spdf)
g.bc = gstat(NULL, "TP", bc.spdf$TP ~ 1, bc.spdf)
## Add an additional component to each object. Now they have both TP and Flow
g.log = gstat(g.log, "FLOW", log.spdf$FLOW ~ 1, log.spdf)
g.bc = gstat(g.bc, "FLOW", bc.spdf$FLOW ~ 1, bc.spdf)
## Calculates auto and cross-semivariograms (for plotting against the fitted model).
v.log = variogram(g.log,cutoff=200,width=200/20)
v.bc = variogram(g.bc,cutoff=200,width=200/20)

## Have a quick look at what we have so far. Thers no real differnce except for the ranges on the y axis.
plot(v.log,g.log)
win.graph()
plot(v.bc,g.bc)

##################################################
### The next few steps are designed to use the ###
### geoR function (likfit) to fit a variogram  ###
### model to each individual variable using    ###
### REML. Slower than gstats, fit.variogram.r- ###
### eml function but I'm tired, and it does t- ###
### he same.                                   ###
##################################################

## Create 4 geodata objects one four each variable and transformation.
qual.tp.log<-as.geodata(data.log,coords.col=1:2,data.col=3)
qual.flow.log<-as.geodata(data.log,coords.col=1:2,data.col=4)
qual.tp.bc<-as.geodata(data.bc,coords.col=1:2,data.col=3)
qual.flow.bc<-as.geodata(data.bc,coords.col=1:2,data.col=4)

## Perform likfit on each different variable. Very time consuming.
log.tp.likfit <- likfit(qual.tp.log,ini=c(0.5,0.5),lik.method="REML") 
log.flow.likfit <- likfit(qual.flow.log,ini=c(0.5,0.5),lik.method="REML")
bc.tp.likfit <- likfit(qual.tp.bc,ini=c(0.5,0.5),lik.method="REML")
bc.flow.likfit <- likfit(qual.flow.bc,ini=c(0.5,0.5),lik.method="REML")
## Rip the psill from the likfits
log.tp.psill <- summary(log.tp.likfit)$spatial.component[[2]][1]
log.flow.psill <- summary(log.flow.likfit)$spatial.component[[2]][1]
bc.tp.psill <- summary(bc.tp.likfit)$spatial.component[[2]][1]
bc.flow.psill <- summary(bc.flow.likfit)$spatial.component[[2]][1]
## Rip the nuggets from the likfits
log.tp.nug <- summary(log.tp.likfit)$nugget.component[[2]]
log.flow.nug <- summary(log.flow.likfit)$nugget.component[[2]]
bc.tp.nug <- summary(bc.tp.likfit)$nugget.component[[2]]
bc.flow.nug <- summary(bc.flow.likfit)$nugget.component[[2]]
## Rip the range values from the likfits
log.tp.range <- summary(log.tp.likfit)$practicalRange
log.flow.range <- summary(log.flow.likfit)$practicalRange
bc.tp.range <- summary(bc.tp.likfit)$practicalRange
bc.flow.range <- summary(bc.flow.likfit)$practicalRange

##################################################
### Finally we can fit the LMCR. There is thr- ###
### ee steps in this process. These three ste- ###
### ps are performed for both transformations. ###
##################################################

##Create gstat objects of the data
log.spdf <- SpatialPointsDataFrame(data.log[,1:2],data.log)
bc.spdf <- SpatialPointsDataFrame(data.bc[,1:2],data.bc)

## Calculates auto and cross-semivariograms (for plotting against the fitted model).
v.log = variogram(g.log,cutoff=200,width=200/20)
v.bc = variogram(g.bc,cutoff=200,width=200/20)

### STEP 1. Create A gstat object with TP and flow.
##The log transformed dataset
g.log = gstat(NULL, "TP", log.spdf$TP ~ 1, log.spdf)
g.log = gstat(g.log, "FLOW", log.spdf$FLOW ~ 1, log.spdf)
##The box.cox transformed dataset
g.bc = gstat(NULL, "TP", bc.spdf$TP ~ 1, bc.spdf)
g.bc = gstat(g.bc, "FLOW", bc.spdf$FLOW ~ 1, bc.spdf)



### STEP 2. Add the three vgm components determined using likfit.
##The logged transformed dataset (fill.all will overwrite the vgm settigs for all variograms)
g.log = gstat(g.log,id="TP",model=vgm(log.tp.psill,"Sph",log.tp.range,log.tp.nug))
g.log = gstat(g.log,id="FLOW",model=vgm(log.flow.psill,"Sph",log.flow.range,log.flow.nug))
##The box.cox transformed dataset (fill.all will overwrite the vgm settigs for all variograms)
g.bc = gstat(g.bc,id="TP",model=vgm(bc.tp.psill,"Sph",bc.tp.range,bc.tp.nug))
g.bc = gstat(g.bc,id="FLOW",model=vgm(bc.flow.psill,"Sph",bc.flow.range,bc.flow.nug))

##################################################
### The next few steps are designed to use the ###
### geoR function (variofit) to fit a variogr- ###
### am model to the cross-variogram. This isn- ###
### 't nessarcy, but makes things easier.      ###
##################################################
## Calcuate the cross-semivariogram clouds for geoR.
v.log.cloud = variogram(g.log,cloud=TRUE)
v.bc.cloud = variogram(g.bc,cloud=TRUE)

## Create a dummy geodata object
dummy.geodata<-as.geodata(data.log,coords.col=1:2,data.col=3)
## Create a dummy variog object using dummy values (one for log and one for flow)
dummy.var1<-variog(qual.tp.log,estimator.type="classical",option=c("cloud"))
dummy.var2<-variog(qual.tp.bc,estimator.type="classical",option=c("cloud"))
## Replace the dummy values with the cross-variogram values
dummy.var1$u<-v.log.cloud$dist
dummy.var1$v<-v.log.cloud$gamma
dummy.var2$u<-v.bc.cloud$dist
dummy.var2$v<-v.bc.cloud$gamma
## Use geoR's variofit to fit a variogram model to the cross-variogram
## (doesn't need to be great) as these are just initial values for the
## vgm component of the final gstat objects.
log<-variofit(dummy.var1,ini=c(0.5,0.5),cov.model="spherical", fix.nugget=F,max.dist=200,minimisation.function="nls")
bc<-variofit(dummy.var2,ini=c(0.5,0.5),cov.model="spherical", fix.nugget=F,max.dist=200,minimisation.function="nls")

##########################################
##########################################
########## BACK TO THE FITTING. ##########
##########################################
##########################################
g.log = gstat(g.log,id=c("TP","FLOW"),model=vgm(summary(log)$spatial.component[[1]][1], "Sph", log$practicalRange, log$nugget))
g.log.fill.all = gstat(g.log,id=c("TP","FLOW"),model=vgm(summary(log)$spatial.component[[1]][1], "Sph", log$practicalRange, log$nugget),fill.all=TRUE)
g.bc = gstat(g.bc,id=c("TP","FLOW"),model=vgm(summary(bc)$spatial.component[[1]][1], "Sph", bc$practicalRange, bc$nugget))
g.bc.fill.all = gstat(g.bc,id=c("TP","FLOW"),model=vgm(summary(bc)$spatial.component[[1]][1], "Sph", bc$practicalRange, bc$nugget),fill.all=TRUE)

### STEP 3. Create the fits using fit.lmc
##The log transformed fits based on all likfit output (fill.all=FALSE)
log.fit=fit.lmc(v.log,g.log)
log.fit1.T = fit.lmc(v.log, g.log,fit.ranges=TRUE,fit.lmc=FALSE,fit.method=1,warn.if.neg=TRUE)
log.fit2.T = fit.lmc(v.log, g.log,fit.ranges=TRUE,fit.lmc=FALSE,fit.method=2,warn.if.neg=TRUE)
log.fit6.T = fit.lmc(v.log, g.log,fit.ranges=TRUE,fit.lmc=FALSE,fit.method=6,warn.if.neg=TRUE)
log.fit7.T = fit.lmc(v.log, g.log,fit.ranges=TRUE,fit.lmc=FALSE,fit.method=7,warn.if.neg=TRUE)
log.fit.fill.all.T=fit.lmc(v.log,g.log.fill.all,fit.ranges=TRUE,fit.lmc=FALSE)
log.fit1.F = fit.lmc(v.log, g.log,fit.ranges=FALSE,fit.lmc=FALSE,fit.method=1,warn.if.neg=TRUE)
log.fit2.F = fit.lmc(v.log, g.log,fit.ranges=FALSE,fit.lmc=FALSE,fit.method=2,warn.if.neg=TRUE)
log.fit6.F = fit.lmc(v.log, g.log,fit.ranges=FALSE,fit.lmc=FALSE,fit.method=6,warn.if.neg=TRUE)
log.fit7.F = fit.lmc(v.log, g.log,fit.ranges=FALSE,fit.lmc=FALSE,fit.method=7,warn.if.neg=TRUE)
log.fit.fill.all.F=fit.lmc(v.log,g.log.fill.all,fit.ranges=FALSE,fit.lmc=FALSE)
##Check the determinates of the model
det.fit(log.fit)
det.fit(log.fit1.T)
det.fit(log.fit2.T)
det.fit(log.fit6.T)
det.fit(log.fit7.T)
det.fit(log.fit.fill.all.T)
det.fit(log.fit1.F)
det.fit(log.fit2.F)
det.fit(log.fit6.F)
det.fit(log.fit7.F)
det.fit(log.fit.fill.all.F)

##Plotting the results. I'm only plotting the one that worked and the following.
plot(v.log,log.fit,main="fit.ranges=FALSE")
win.graph()
plot(v.log,log.fit7.T,main="fit.ranges=TRUE")

##The box.cox transformed fits.
bc.fit=fit.lmc(v.bc,g.bc)
bc.fit1.T = fit.lmc(v.bc, g.bc,fit.ranges=TRUE,fit.lmc=FALSE,fit.method=1,warn.if.neg=TRUE)
bc.fit2.T = fit.lmc(v.bc, g.bc,fit.ranges=TRUE,fit.lmc=FALSE,fit.method=2,warn.if.neg=TRUE)
bc.fit6.T = fit.lmc(v.bc, g.bc,fit.ranges=TRUE,fit.lmc=FALSE,fit.method=6,warn.if.neg=TRUE)
bc.fit7.T = fit.lmc(v.bc, g.bc,fit.ranges=TRUE,fit.lmc=FALSE,fit.method=7,warn.if.neg=TRUE)
bc.fit.fill.all.T=fit.lmc(v.bc,g.bc.fill.all,fit.ranges=TRUE,fit.lmc=FALSE)
bc.fit1.F = fit.lmc(v.bc, g.bc,fit.ranges=FALSE,fit.lmc=FALSE,fit.method=1,warn.if.neg=TRUE)
bc.fit2.F = fit.lmc(v.bc, g.bc,fit.ranges=FALSE,fit.lmc=FALSE,fit.method=2,warn.if.neg=TRUE)
bc.fit6.F = fit.lmc(v.bc, g.bc,fit.ranges=FALSE,fit.lmc=FALSE,fit.method=6,warn.if.neg=TRUE)
bc.fit7.F = fit.lmc(v.bc, g.bc,fit.ranges=FALSE,fit.lmc=FALSE,fit.method=7,warn.if.neg=TRUE)
bc.fit.fill.all.F=fit.lmc(v.bc,g.bc.fill.all,fit.ranges=FALSE,fit.lmc=FALSE)
##Check the determinates of the model
det.fit(bc.fit)
det.fit(bc.fit1.T)
det.fit(bc.fit2.T)
det.fit(bc.fit6.T)
det.fit(bc.fit7.T)
det.fit(bc.fit.fill.all.T)
det.fit(bc.fit1.F)
det.fit(bc.fit2.F)
det.fit(bc.fit6.F)
det.fit(bc.fit7.F)
det.fit(bc.fit.fill.all.F)

##Plotting the results. I'm only plotting the one that worked and the following.

plot(v.bc,bc.fit,main="fit.ranges=FALSE")
win.graph()
plot(v.bc,bc.fit7.T,main="fit.ranges=TRUE")

##Compare the values for the four best models.
log.fit
log.fit7.T
bc.fit
bc.fit7.T
