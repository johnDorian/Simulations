#!/usr/bin/Rscript
################################################################################
################################################################################
### Title: LMCR.R ##############################################################
### Author: Jason Lessels j.lessels@usyd.edu.au ################################
### Created: Tue Jan 26 16:39:42 2010 ##########################################
### Modified: Tue Jan 26 16:39:42 2010 #########################################
### Aim: To fit a one dimensional LMCR #########################################
### Description: This will use the gstat package to fit a LMCR #################
################################################################################
### NOTES: This will use ML not REML ###########################################
################################################################################
### TODO: Transfer the script from it's source to be compable with this code. ##
################################################################################

## Open the directory for the data

setwd("~/Documents/code/Simulations/")

data<-read.csv("TP_flow.csv", header=TRUE)
## backup the data
back.up <- data
## fix up the data
data <- data[,-c(1,2)]
names(data) <- c("X","Y","TP","FLOW")
data$TP <- log(data$TP+1)
data$FLOW <- log(data$FLOW+1)
##Non-spatial statistics


library(gstat)

##Creates spatial points class

qual<-SpatialPointsDataFrame(data[,1:2],data)

summary(qual)

plot(qual)

###MODELLING AUTO-SEMIVARIOGRAMS

#Calculates experimental auto-semivariogram

b0.vgm <- variogram(qual$TP ~ 1, qual)
b1.vgm <- variogram(qual$FLOW ~ 1, qual)


quartz()
par(mfrow=c(2,1))

plot(b0.vgm[,2],b0.vgm[,3])
plot(b1.vgm[,2],b1.vgm[,3])

#Fits auto-semivariograms

b0.fit <- fit.variogram(b0.vgm,model=vgm(0,"Sph",30000,10))
b1.fit <- fit.variogram(b1.vgm,model=vgm(0,"Sph",30000,10))


quartz()
par(mfrow=c(2,1))

plot(b0.vgm,b0.fit)
plot(b1.vgm,b1.fit)

###MODELLING LINEAR MODEL OF COERGIONALISATION

#Creates a dataset for auto-semivariogram estimation
g = gstat(NULL, "TP", qual$TP ~ 1, qual)
#Somehow links b1 with b0 so it knows cross-semivariograms are needed
g = gstat(g, "FLOW", qual$FLOW ~ 1, qual)
#Calculates auto and cross-semivariograms
v = variogram(g,cutoff=200,width=200/20)
plot(v)
#Defines form of LMCR - range, variogram model
g = gstat(g, model = vgm(1, "Sph", 30000, 1), fill.all = TRUE)
#Fits LMCR model
g.fit = fit.lmc(v, g)

plot(v, g.fit)







