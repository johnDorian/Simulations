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






##################################################
### Finally we can fit the LMCR. There is thr- ###
### ee steps in this process. These three ste- ###
### ps are performed for both transformations. ###
##################################################
##Create gstat objects of the data
log.spdf <- SpatialPointsDataFrame(data.log[,1:2],data.log)
bc.spdf <- SpatialPointsDataFrame(data.bc[,1:2],data.bc)

### STEP 1. Create A gstat object with TP and flow.
##The log transformed dataset
g.log = gstat(NULL, "TP", log.spdf$TP ~ 1, log.spdf,,maxdist=200)
g.log = gstat(g.log, "FLOW", log.spdf$FLOW ~ 1, log.spdf,maxdist=200)
##The box.cox transformed dataset
g.bc = gstat(NULL, "TP", bc.spdf$TP ~ 1, bc.spdf,,maxdist=200)
g.bc = gstat(g.bc, "FLOW", bc.spdf$FLOW ~ 1, bc.spdf,maxdist=200)
#Use likfit to estimate the variable of the variogram models. (This has been done and is in earlier versions of this code, but the data file is quicker to use.)
load("likfit.Rdata")
#variogram models
log.tp.vgm=vgm(log.tp.psill,"Mat",log.tp.range,log.tp.nug)
log.flow.vgm=vgm(log.flow.psill,"Mat",log.flow.range,log.flow.nug)
bc.tp.vgm=vgm(bc.tp.psill,"Mat",bc.tp.range,bc.tp.nug)
bc.flow.vgm=vgm(bc.flow.psill,"Mat",bc.flow.range,bc.flow.nug)
##The logged transformed dataset (fill.all will overwrite the vgm settigs for all variograms anyway but just incase)
g.log = gstat(g.log,id="TP",model=log.tp.vgm,maxdist=200)
g.log = gstat(g.log,id="FLOW",model=log.flow.vgm,maxdist=200)
##The box.cox transformed dataset (fill.all will overwrite the vgm settigs for all variograms)
g.bc = gstat(g.bc,id="TP",model=bc.tp.vgm,maxdist=200)
g.bc = gstat(g.bc,id="FLOW",model=bc.flow.vgm,maxdist=200)
## Calculates auto and cross-semivariograms (for plotting against the fitted model).
v.log = variogram(g.log,cutoff=200,width=200/20)
v.bc = variogram(g.bc,cutoff=200,width=200/20)
##################################################
### The next few steps are designed to use the ###
### geoR function (variofit) to fit a variogr- ###
### am model to the cross-variogram. This isn- ###
### 't nessarcy, but makes things easier.      ###
##################################################
## Create 4 geodata objects one four each variable and transformation.
qual.tp.log<-as.geodata(data.log,coords.col=1:2,data.col=3)
qual.flow.log<-as.geodata(data.log,coords.col=1:2,data.col=4)
qual.tp.bc<-as.geodata(data.bc,coords.col=1:2,data.col=3)
qual.flow.bc<-as.geodata(data.bc,coords.col=1:2,data.col=4)
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
log<-variofit(dummy.var1,ini=c(0.5,0.5),cov.model="matern",fix.nugget=F,minimisation.function="nls",weights="equal")
bc<-variofit(dummy.var2,ini=c(0.5,27),cov.model="matern", fix.nugget=F,minimisation.function="nls",weights="equal")
##############################

###Create a vgm using the variogram cloud values.
log.vgm <- vgm(summary(log)$spatial.component[[1]][1], "Mat", log$practicalRange, log$nugget)
bc.vgm <- vgm(summary(bc)$spatial.component[[1]][1], "Mat", bc$practicalRange, bc$nugget)
#The cross covariance models
g.log = gstat(g.log,id=c("TP","FLOW"),maxdist=200,model=log.vgm,fill.all=TRUE)
g.bc = gstat(g.bc,id=c("TP","FLOW"),model=bc.vgm,fill.all=TRUE)

### STEP 3. Create the fits using fit.lmc
##The log transformed fits based on all likfit output (fill.all=TRUE)
log.fit=fit.lmc(v.log,g.log,fit.lmc=TRUE,fit.ranges=FALSE)
bc.fit=fit.lmc(v.bc,g.bc,fit.lmc=TRUE,fit.ranges=FALSE)
##Check the determinates of the model
det.fit(bc.fit)
det.fit(log.fit)
##Plotting the results.
plot(v.bc,bc.fit,main="boxcox transformation with matern cov model")
win.graph()
plot(v.bc,bc.fit,main="log transformation with matern cov model")
##Compare the values for the four best models.
log.fit
bc.fit
