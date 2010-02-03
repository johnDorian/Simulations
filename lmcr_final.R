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
## Now so we can easly check the determinates
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
### The next few steps are designed to use the ###
### geoR function (variofit) to fit a variogr- ###
### am model to the cross-variogram. This isn- ###
### 't nessarcy, but makes things easier.      ###
##################################################
## Calcuate the cross-semivariogram clouds for geoR.
v.log.cloud = variogram(g.log,cloud=TRUE)
v.bc.cloud = variogram(g.log,cloud=TRUE)

## Create a dummy geodata object
dummy.geodata<-as.geodata(data.log,coords.col=1:2,data.col=3)
## Create a dummy variog object using dummy values (one for log and one for flow)
dummy.var1<-variog(qual.tp.log,estimator.type="classical",option=c("cloud"))
dummy.var2<-variog(qual.tp.log,estimator.type="classical",option=c("cloud"))
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


#######################################


### STEP 2. Add the three vgm components determined using likfit.
##The logged transformed dataset (fill.all will overwrite the vgm settigs for all variograms)
g.log = gstat(g.log,id="TP",model=vgm(log.tp.psill,"Sph",log.tp.range,log.tp.nug))
g.log = gstat(g.log,id="FLOW",model=vgm(log.flow.psill,"Sph",log.flow.range,log.flow.nug))
g.log = gstat(g.log,id=c("TP","FLOW"),model=vgm(summary(log)$spatial.component[[1]][1], "Sph", log$practicalRange, log$nugget))
g.log.fill.all = gstat(g.log,id=c("TP","FLOW"),model=vgm(summary(log)$spatial.component[[1]][1], "Sph", log$practicalRange, log$nugget),fill.all=TRUE)
##The box.cox transformed dataset (fill.all will overwrite the vgm settigs for all variograms)
g.bc = gstat(g.bc,id="TP",model=vgm(bc.tp.psill,"Sph",bc.tp.range,bc.tp.nug))
g.bc = gstat(g.bc,id="FLOW",model=vgm(bc.flow.psill,"Sph",bc.flow.range,bc.flow.nug))
g.bc = gstat(g.bc,id=c("TP","FLOW"),model=vgm(summary(bc)$spatial.component[[1]][1], "Sph", bc$practicalRange, bc$nugget))
g.bc.fill.all = gstat(g.bc,id=c("TP","FLOW"),model=vgm(summary(bc)$spatial.component[[1]][1], "Sph", bc$practicalRange, bc$nugget),fill.all=TRUE)

### STEP 3. Create the fits using fit.lmc
##The log transformed fits.
log.fit1 = fit.lmc(v.log, g.log,fit.sills=TRUE,fit.ranges=TRUE,fit.lmc=TRUE,fit.method=1,warn.if.neg=TRUE)
log.fit2 = fit.lmc(v.log, g.log,fit.sills=TRUE,fit.ranges=TRUE,fit.lmc=TRUE,fit.method=2,warn.if.neg=TRUE)
log.fit6 = fit.lmc(v.log, g.log,fit.sills=TRUE,fit.ranges=TRUE,fit.lmc=TRUE,fit.method=6,warn.if.neg=TRUE)
log.fit7 = fit.lmc(v.log, g.log,fit.sills=TRUE,fit.ranges=TRUE,fit.lmc=TRUE,fit.method=7,warn.if.neg=TRUE)
log.fit=fit.lmc(v.log,g.log,fit.ranges=TRUE,fit.lmc=FALSE)
log.fit.fill.all=fit.lmc(v.log,g.log.fill.all,fit.ranges=TRUE,fit.lmc=FALSE)
det.fit(log.fit)
det.fit(log.fit.fill.all)
plot(v.log,log.fit.fill.all,main="fill all")
win.graph()
plot(v.log,log.fit,main="using likfit")

##Check the determinate using function (at start).
det.fit(log.fit1)
det.fit(log.fit2)
det.fit(log.fit6)
det.fit(log.fit7)
plot(v.log, log.fit1)
plot(v.log, log.fit2)
plot(v.log, log.fit6)
plot(v.log, log.fit7)
log.fit7

##The box.cox transformed fits.
bc.fit1 = fit.lmc(v.bc, g.bc,fit.lmc=FALSE,fit.method=1)
bc.fit2 = fit.lmc(v.bc, g.bc,fit.lmc=TRUE,fit.method=2)
bc.fit6 = fit.lmc(v.bc, g.bc,fit.lmc=TRUE,fit.method=6)
bc.fit7 = fit.lmc(v.bc, g.bc,fit.ranges=TRUE) #The default fit.method is option 7.
bc.fit7.2 = fit.lmc(v.bc, g.bc)
##Check the determinate using function (at start).
det.fit(bc.fit1)
det.fit(bc.fit2)
det.fit(bc.fit6)
det.fit(bc.fit7)
plot(v.bc, bc.fit1)
plot(v.bc, bc.fit2)
plot(v.bc, bc.fit6)
plot(v.bc, bc.fit7)
win.graph()
plot(v.bc,bc.fit7.2)











## Now using the values fitted using variofit set the vgm() in the gstat() function for the cross-variogram.
g.log = gstat(g.log,id=c("TP","FLOW"), model = vgm(summary(log)$spatial.component[[1]][1], "Sph", log$practicalRange, log$nugget), fill.all = TRUE)
g.bc = gstat(g.bc, model = vgm(summary(bc)$spatial.component[[1]][1], "Sph", bc$practicalRange, bc$nugget), fill.all = TRUE)

## Using the gstat objects, and the variograms use fit.lmc to fit a final model to the variograms.
## IMPORTANT: READ THE DOCUMENTATION. As defualt fit.ranges=FALSE,
## and fit.lmc=!fit.ranges. If fit.lmc=TRUE, than each partial si-
## ll will be  a positive definite.

## A function to test the determenant of each fit. Only works for this
## data because of the names within the function (i.e. TP and FLOW).
###Now to test if the deteriment equals zero.
det.fit <- function(object){
tpn <- object$model$TP$psill[1]
fln <- object$model$FLOW$psill[1]
tpfln <- object$model$TP.FLOW$psill[1]
d <- matrix(c(tpn,tpfln,tpfln,fln),c(2,2))
print(det(d))
}
## Now the fit.lmc models.
log.fit.1 = fit.lmc(v.log, g.log,fit.sills=TRUE,fit.ranges=TRUE,fit.lmc=FALSE,fit.method=1)
log.fit.2 = fit.lmc(v.log, g.log,fit.sills=TRUE,fit.ranges=TRUE,fit.lmc=FALSE,,fit.method=2)
log.fit.6 = fit.lmc(v.log, g.log,fit.sills=TRUE,fit.ranges=TRUE,fit.lmc=FALSE,fit.method=6)
log.fit.7 = fit.lmc(v.log, g.log,fit.sills=TRUE,fit.ranges=TRUE,fit.lmc=FALSE,,fit.method=7)
det.fit(log.fit.1)
#0.2850109
det.fit(log.fit.2)
#0.1191131
det.fit(log.fit.6)
#0.3860414
det.fit(log.fit.7)
#0.08932527
bc.fit1 = fit.lmc(v.bc, g.bc,fit.sills=TRUE,fit.ranges=TRUE,fit.lmc=FALSE,fit.method=1)
bc.fit2 = fit.lmc(v.bc, g.bc,fit.sills=TRUE,fit.ranges=TRUE,fit.lmc=FALSE,fit.method=2)
bc.fit6 = fit.lmc(v.bc, g.bc,fit.sills=TRUE,fit.ranges=TRUE,fit.lmc=FALSE,fit.method=6)
bc.fit7 = fit.lmc(v.bc, g.bc,fit.sills=TRUE,fit.ranges=TRUE,fit.lmc=FALSE,fit.method=7)
det.fit(bc.fit1)
# 0.0978011
det.fit(bc.fit2)
# 0.07568633
det.fit(bc.fit6)
# 0.1263060
det.fit(bc.fit7)
# 0.05446267



