####Modelling coregionalised data####

#Uses GSAT
##Start:  6th April 2009
##Finish: 6th April 2000


setwd("C:/Tom/Sync/Research/Modelling/Rainfall/Data")

soil<-read.table("rainfall.txt", header=TRUE, sep="\t", dec=".")

##Non-spatial statistics

hist(soil$b0)

library(gstat)

##Creates spatial points class

cash<-SpatialPointsDataFrame(soil[,1:2],soil)

summary(cash)

plot(cash)

###MODELLING AUTO-SEMIVARIOGRAMS

#Calculates experimental auto-semivariogram

b0.vgm <- variogram(cash$b0 ~ 1, cash)
b1.vgm <- variogram(cash$b1 ~ 1, cash)
resb0.vgm <- variogram(cash$resb0 ~ 1, cash)
resb1.vgm <- variogram(cash$resb1 ~ 1, cash)


win.graph()
par(mfrow=c(2,2))

plot(b0.vgm[,2],b0.vgm[,3])
plot(b1.vgm[,2],b1.vgm[,3])
plot(resb0.vgm[,2],resb0.vgm[,3])
plot(resb1.vgm[,2],resb1.vgm[,3])

#Fits auto-semivariograms

b0.fit <- fit.variogram(b0.vgm,model=vgm(0,"Sph",30000,10))
b1.fit <- fit.variogram(b1.vgm,model=vgm(0,"Sph",30000,10))
resb0.fit <- fit.variogram(resb0.vgm,model=vgm(0,"Sph",30000,10))
resb1.fit <- fit.variogram(resb1.vgm,model=vgm(0,"Sph",30000,10))

win.graph()
par(mfrow=c(2,2))

plot(b0.vgm,b0.fit)
plot(b1.vgm,b1.fit)
plot(resb0.vgm,resb0.fit)
plot(resb1.vgm,resb1.fit)

###MODELLING LINEAR MODEL OF COERGIONALISATION

#Creates a dataset for auto-semivariogram estimation
g = gstat(NULL, "b0", cash$b0 ~ 1, cash)
#Somehow links b1 with b0 so it knows cross-semivariograms are needed
g = gstat(g, "b1", cash$b1 ~ 1, cash)
#Calculates auto and cross-semivariograms
v = variogram(g)

#Defines form of LMCR - range, variogram model
g = gstat(g, model = vgm(1, "Sph", 30000, 1), fill.all = TRUE)
#Fits LMCR model
g.fit = fit.lmc(v, g)

plot(v, g.fit)



