################################################################################
################################################################################
### Title: LMCR_FINAL.R ########################################################
### Author: Jason Lessels j.lessels@usyd.edu.au ################################
### Created: Wed Feb  22/03/2010 ###############################################
### Modified: Wed Feb  30/03/2010 ##############################################
### Aim: To perform unconditional simulation using an lmcr model ###############
### Description:This script performs a lmcr using two main steps, the first step
### is by using the gstat library to get the cross co variogram and the second
### step is by using the method of Lark 2003 to use simulated annealing process
### to fit the model. Than using the gstat package unconditional simulations are### conducted.
################################################################################
### NOTES: This script performs unconditional simulation of water quality and discharge data.
################################################################################
### TODO:  Once everthing is working to merge this with the lmcr.R script.######
################################################################################

### Set the working directory
setwd("~/Documents/code/Simulations/")
### Load the datasets.
load("tp_flow_krige.Rdata")
### Load the gstat library
library(gstat);library(TeachingDemos);library(geoR)
### Load the lmcr functions - from script.
source("lmcr_convert.R")
## Clean up the dataset
names(tp_flow)
data <-tp_flow[,-1]
head(data)
data <- subset(data,!is.nan(data)&!is.na(data))
data <- data.frame(X=data[,3],Y=1,TP=data[,2],FLOW=data[,1])
head(data)
flow.lambda <- boxcox.fit(data$FLOW+0.01)$lambda
tp.lambda <- boxcox.fit(data$TP+0.001)$lambda
### Transform the data using a log scale transformation.
data$TP <- bct(data$TP+0.001,tp.lambda)
data$FLOW <- bct(data$FLOW+0.01,flow.lambda)

### Create gstat object of the data
spdf <- SpatialPointsDataFrame(data[,1:2],data)
### Create A gstat object with TP and flow.
g = gstat(NULL, "TP", TP ~ 1,spdf,maxdist=200,dummy=TRUE,beta=mean(data$TP))
g = gstat(g, "FLOW", FLOW ~ 1,spdf,maxdist=200,dummy=TRUE,beta=mean(data$FLOW))


### Create cross and auto variograms of the gstat object, using width to set the bins.
v = variogram(g,cutoff=200,width=200/20)

### Plot the cross covariogram and get some initial values for it.
plot(v)
############################################################
#### Set and get the parameters for the lmcr function. #####
############################################################

### Set the temperature of the annealing process. See Lark 2003 for more details.
cpar=30000
### Set the model type that is desired. See the gammah function for a translation. ONLY 4 HAS BEEN TESTED (exp). The others should work.
modtyp=4 
### Get the gstat estimated values of the variogram structure to initialise the lmcr.
### Order from gstat object is second variable, cross, cross, first variable.
covar<-c(30,0.05,0.002)
### Now for the guesses. (You have to give the distance twice as the function also handles a double exp model). This dosent effect the results
guessa<-50
### Run the lmcr function
model<-lmcr(g,v,covar,guessa,modtyp,cpar)
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

#####Simulation zeit###
####Get the data organised...
###Make a grid for simulation. (In this example I want an hourly grid using days as the units for the x column and a constant 1 for the y column.
x=(seq(1,700)-1)*0.04166667
###Run the simulation using the gstat model from the lmcr function using the simulation grid to simulate to and nsim=no of simulations to perform.
f<-GaussRF(x=x, model="exp",param=c(mean(data$TP),9.533031e-04 ,3.987555e-06),n=1000)

###Now see if all this has worked
par(mfrow=c(2,1))
plot(sim1$FLOW.sim1,type="l",ylab="FLOW")
plot(sim1$TP.sim1,type="l",ylab="TP")### 

