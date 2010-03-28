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

### Set the working directory
setwd("~/Documents/code/Simulations/")
### Load the datasets.
data<-read.csv("TP_flow.csv", header=TRUE)
### Load the gstat library
library(gstat);library(TeachingDemos);library(geoR)
### Load the lmcr functions - from script.
source("lmcr_convert.R")
## Clean up the dataset
data <- data[,-c(1,2)]
names(data) <- c("X","Y","TP","FLOW")
### Transform the data using a log scale transformation.
data$TP <- bct(data$TP+0.001,boxcox.fit(data$TP+0.001)$lambda)
data$FLOW <- bct(data$FLOW+1,boxcox.fit(data$FLOW+1)$lambda)

### Create gstat object of the data
spdf <- SpatialPointsDataFrame(data[,1:2],data)

### Create A gstat object with TP and flow.
g = gstat(NULL, "TP", spdf$TP ~ 1,spdf,maxdist=200)
g = gstat(g, "FLOW", spdf$FLOW ~ 1,spdf,maxdist=200)


### Create cross and auto variograms of the gstat object, using width to set the bins.
v = variogram(g,cutoff=200,width=200/20)

### Plot the cross covariogram and get some initial values for it.
plot(v)
############################################################
#### Set and get the parameters for the lmcr function. #####
############################################################

### If distance should be the same. 1=guessa is correct. 0= fit range.
lock=0
### Set the temperature of the annealing process. See Lark 2003 for more details.
cpar=30000
### Set the model type that is desired. See the gammah function for a translation. ONLY 4 HAS BEEN TESTED (exp). The others should work.
modtyp=4 
### Set the constraints of the returned matrix to be positive definitive.
icvp=1
### Set the weighting option. See the fcn function for more details. wgt=1 gives equal weighting for the weighted sums-of-squares
wgt=1
### Set the maximum distance of the variogram - based on the object v
maxdist=200
### Get the gstat estimated values of the variogram structure to initialise the lmcr.
### Order from gstat object is second variable, cross, cross, first variable.
covar<-matrix(c(500,0.3,0.3,0.002),c(2,2))
### Now for the guesses. (You have to give the distance twice as the function also handles a double exp model). This dosent effect the results
guessa<-c(50,50)
### Run the lmcr function
model<-lmcr(g,v,wgt,icvp,cpar,modtyp,covar,maxdist,c(50,50),lock)

### Check out the results of the fit. 
plot(model$variogram,model$g)
### Now check the positive defitive of the results (i.e that the det is >0)
det(matrix(c(model$c[4],model$c[5],model$c[5],model$c[6]),c(2,2)))
#1.534631
### Now check the Structural correlation <1 ~0.7 is best.
model$c[5]/(sqrt(model$c[4]*model$c[6]))
#0.5198716
### Save the results. (Just so I don't need to do all this again).
save(file="lmcr.Rdata",model)

#####For the simulations
xy <- expand.grid(1:100, 1)
names(xy) <- c("x","y")
g.dummy <- gstat(formula = z~1, locations = ~x+y, dummy = TRUE, beta = 0,model = vgm(1,"Exp",15), nmax = 20)
yy <- predict(g.dummy, newdata = xy, nsim = 4)
####Get the data organised...
setwd("~
