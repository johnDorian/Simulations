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
library(gstat)
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

### If distance should be the same. 1=leave range alone.
lock=0
### Set the temperature of the annealing process. See Lark 2003 for more details. 
cpar=1500 
### Set the model type that is desired. See the gammah function for a translation. ONLY 4 HAS BEEN TESTED (exp).
modtyp=4 
### Set the constraints of the returned matrix to be positive definitive.
icvp=1
### Set the weighting option. See the fcn function for more details. wgt=1 gives equal weighting for the weighted sums-of-squares
wgt=1
### Set the maximum distance of the variogram - based on the object v
maxdist=200
### Get the gstat estimated values of the variogram structure to initialise the lmcr.
### Order is FLOW,cross,cross,TP
covar<-matrix(c(1,0.3,0.3,0.2),c(2,2))
### make sure covar is the same as semvar
## 	a.1	c
##	c	a.2
### And now we run the model.
test<-lmcr(v,wgt,icvp,cpar,modtyp,covar,maxdist,c(29,29),lock)

### Use the results once everything is done to put back into the 
### The results from the lmcr are in the format (within test$c)
### 1		2		3	4	5	6
### a.1_nugget,c.1_nugget,a.2_nugget,a.1_sill,c.1_sill,a.2_sill
g = gstat(g,id=c("TP","FLOW"),model=vgm(test$c[5],"Exp",test$distance,test$c[2]),locations=~x+y)
g = gstat(g,"TP",model=vgm(test$c[6],"Exp",test$distance,test$c[3]))
g = gstat(g,"FLOW",model=vgm(test$c[4],"Exp",test$distance,test$c[1]))
### Check out the results of the fit. 
plot(v,g)
### Now check the positive defitive of the results (i.e that the det is >0)
det(matrix(c(test$c[4],test$c[5],test$c[5],test$c[6]),c(2,2)))
#1.534631
### Now check the Structural correlation <1 ~0.7 is best.
test$c[5]/(sqrt(test$c[4]*test$c[6]))
#0.5198716


#####For the simulations
xy <- expand.grid(1:100, 1)
names(xy) <- c("x","y")
g.dummy <- gstat(formula = z~1, locations = ~x+y, dummy = TRUE, beta = 0,model = vgm(1,"Exp",15), nmax = 20)
yy <- predict(g.dummy, newdata = xy, nsim = 4)
