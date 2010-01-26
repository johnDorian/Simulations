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

setwd("C:/Tom/Sync/Research/Modelling/Rainfall/Data")

data<-read.table("rainfall.txt", header=TRUE, sep="\t", dec=".")

##Non-spatial statistics

hist(data$)

library(gstat)

##Creates spatial points class

qual<-SpatialPointsDataFrame(data[,1:2],data)

summary(qual)

plot(qual)

###MODELLING AUTO-SEMIVARIOGRAMS

#Calculates experimental auto-semivariogram

b0.vgm <- variogram(qual$b0 ~ 1, cash)
b1.vgm <- variogram(qual$b1 ~ 1, cash)


quartz()
par(mfrow=c(2,2))

plot(b0.vgm[,2],b0.vgm[,3])
plot(b1.vgm[,2],b1.vgm[,3])

#Fits auto-semivariograms

b0.fit <- fit.variogram(b0.vgm,model=vgm(0,"Sph",30000,10))
b1.fit <- fit.variogram(b1.vgm,model=vgm(0,"Sph",30000,10))


quartz()
par(mfrow=c(2,2))

plot(b0.vgm,b0.fit)
plot(b1.vgm,b1.fit)

###MODELLING LINEAR MODEL OF COERGIONALISATION

#Creates a dataset for auto-semivariogram estimation
g = gstat(NULL, "b0", qual$b0 ~ 1, cash)
#Somehow links b1 with b0 so it knows cross-semivariograms are needed
g = gstat(g, "b1", qual$b1 ~ 1, cash)
#Calculates auto and cross-semivariograms
v = variogram(g)

#Defines form of LMCR - range, variogram model
g = gstat(g, model = vgm(1, "Sph", 30000, 1), fill.all = TRUE)
#Fits LMCR model
g.fit = fit.lmc(v, g)

plot(v, g.fit)







