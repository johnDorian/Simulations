################################################################################
################################################################################
### Title: LMCR_FINAL.R ########################################################
### Author: Jason Lessels j.lessels@usyd.edu.au ################################
### Created: Wed Feb  22/30/2010 ###############################################
### Modified: Wed Feb  22/03/2010 ##############################################
### Aim: To perform LMCR #######################################################
### Description:This script performs a lmcr using two main steps, the first step
### is by using the gstat library to get the cross co variogram and the second
### step is by using the method of Lark 2003 to use simulated annealing process
### to fit the model.
################################################################################
### NOTES: This script uses the lmcr_covert.R script for the Larks methods. The
### main purpose of this script is to provide a cleaner easier to understand 
### method of the lmcr.R script.
################################################################################
### TODO:  Once everthing is working to merge this with the lmcr.R script.######
################################################################################

###GET SOME USEFUL FUNCTIONS OUT OF THE WAY.
### Let's universalise plotting on all Oprating Systems.
win.graph = function(){if(.Platform$OS=='unix') x11() else win.graph=win.graph}
### Set the working directory
setwd("~/Documents/code/Simulations/")
### Load the datasets.
data<-read.csv("TP_flow.csv", header=TRUE)
### Load the nessacry libraries. A couple for variograms, and some for the boxcox transformations.
library(gstat)
### Load the lmcr functions - from script.
source("lmcr_convert.R")
## Clean up the dataset - create a back.up of the data
back.up <- data
data <- data[,-c(1,2)]
names(data) <- c("X","Y","TP","FLOW")
### Transform the data using a log scale transformation.
data$TP <- log(data$TP+0.01)
data$FLOW <- log(data$FLOW+0.01)

### Create gstat objects of the data
spdf <- SpatialPointsDataFrame(data[,1:2],data)

### Create A gstat object with TP and flow.
g = gstat(NULL, "TP", spdf$TP ~ 1, spdf,,maxdist=200)
g = gstat(g, "FLOW", spdf$FLOW ~ 1, spdf,maxdist=200)

### Create variogarms.
v = variogram(g,cutoff=200,width=200/20)

### Plot the cross covariogram and get some initial values for it.
plot(v,g)
#The cross covariance models using initial values, these will change.
g = gstat(g,id=c("TP","FLOW"),model=vgm(0.5,"Exp",48,0.2),fill.all=TRUE)


### Use fit.lmc() function to find the ranges for each covariogram, using gstat to optimise the range values.
fit=fit.lmc(v,g,fit.lmc=FALSE,fit.range=TRUE)

### Set and get the parameters for the lmcr function.

### Get the distance parameters from the gstat function
guessa=ranger(fit)
### Get the amount of observations of each variogram
nlags=matrix(c(length(semvar[,1]),length(semvar[,1]),length(semvar[,1])))
### If distance should be the same. 1=same.
lock=1
### Set the temperature of the annealing process. See Lark 2003 for more details
cpar=200 
### Set the model type that is desired. See the gammah function for a translation
modtyp=4 
### Set the constraints of the returned matrix to be positive definitive.
icvp=1
### Set the weighting option. See the fcn function for more details
wgt=1
### Set the maximum distance of the variogram
maxdist=200
### Get the values of the variogram - bins,gamma etc. Obtain this from the gstat variogram.
semvar<-semVar(v)
### Get the gstat estimated values of the variogram structure to initialise the lmcr.
covar<-covStructure(fit)
### And now we run the model
test<-lmcr(semvar,nlags,wgt,icvp,cpar,modtyp,covar,maxdist,guessa,lock)
