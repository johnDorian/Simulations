################################################################################
################################################################################
### Title: LMCR_FINAL.R ########################################################
### Author: Jason Lessels jason.lessels@sydney.edu.au ##########################
### Created: Wed Feb  22/03/2010 ###############################################
### Modified: Wed Feb  24/04/2010 ##############################################
### Aim: To perform unconditional simulation using an lmcr model ###############
### Description:This script performs a lmcr using two main steps, the first step
### is by using the gstat library to get the cross co variogram and the second
### step is by using the method of Lark 2003 to use simulated annealing process
### to fit the model. The simulation step is performed using a combination of 
### the technique outlined by Bishop/LArk 2006 and the grf function in the geoR 
### package.
################################################################################
### TODO:  Once everthing is working to merge this with the lmcr.R script.######
################################################################################

### Load the datasets.
load("~/Documents/code/Simulations/data/raw_data/tp_flow_krige.Rdata")
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
flow.lambda <- boxcox.fit(data$FLOW+0.001)$lambda
tp.lambda <- boxcox.fit(data$TP+0.001)$lambda
### Transform the data the boxcox formula.
data$TP <- ((data$TP+0.001)^tp.lambda-1)/tp.lambda
data$FLOW <- ((data$FLOW+0.001)^flow.lambda-1)/flow.lambda
### Create gstat object of the data
spdf <- SpatialPointsDataFrame(data[,1:2],data)
### Create A gstat object with TP and flow.
g = gstat(NULL, "Total Phosphorus", TP ~ 1,spdf,maxdist=200,dummy=TRUE,beta=mean(data$TP))
g = gstat(g, "Discharge", FLOW ~ 1,spdf,maxdist=200,dummy=TRUE,beta=mean(data$FLOW))


### Create cross and auto variograms of the gstat object, using width to set the bins.
v = variogram(g,cutoff=150,width=10)

### Plot the cross covariogram and get some initial values for it.
plot(v)
############################################################
#### Set and get the parameters for the lmcr function. #####
############################################################

### Set the temperature of the annealing process. See Lark 2003 for more details.
cpar=20
### Set the model type that is desired. See the gammah function for a translation. ONLY 4 HAS BEEN TESTED (exp). The others should work.
modtyp=4 
### Get the gstat estimated values of the variogram structure to initialise the lmcr. (flow,cross,tp)
covar<-c(5,1.5,0.8)
### Now for the guesse of range (phi) - If a double type model is used, two distances are required.
guessa<-50/3
### Run the lmcr function
model<-lmcr(g,v,covar,guessa,modtyp,cpar)
### Check out the results of the fit. 
plot(model$variogram,model$g)
### Now check the positive defitive of the results (i.e that the det is >0)
det(matrix(c(model$c[4],model$c[5],model$c[5],model$c[6]),c(2,2)))
#2.838906
### Now check the Structural correlation <1 ~0.7 is best.
model$c[5]/(sqrt(model$c[4]*model$c[6]))
# 0.5749334
### Save the results. (Just so I don't need to do all this again).
save(file="lmcr.Rdata",model)
load("lmcr.Rdata")


#####Simulation zeit###
####Get the data organised...
###Make a grid for simulation. (In this example I want an hourly grid using days as the units for the x column and a constant 1 for the y column.
sim=24*365*20
no.sim=2500
x=(seq(1,sim)-1)*0.04166667
grid=cbind(x,1)
###Run the simulation using the gstat model from the lmcr function using the simulation grid to simulate to and nsim=no of simulations to perform.
##cov.pars(psill,range)


###Get the random fields of the sills. - each takes ~ 35 mins for 500, but there was a problem with memory.model$g$model$TP$psill[[2]]
  
###Find the a values for the sill and the nugget using the function 'fiting'
###for the sill - with the order of results a11 = flow, a12 = cross, a21 = cross and a22= tp
source('fiting.R')
a1<-fiting(model$gstat)
###And now the nugget
a0<-fiting(model$gstat,parm="nugget")

#### - start from here (temp. comment for 15/07/10)
sill.1<-grf(sim,grid=grid,nsim=no.sim,cov.pars=c(1,model$gstat$model[[1]][2,3]),nug=0,mean=0)$data
save(sill.1,file="sill1.Rdata")
rm(sill.1)
gc()

sill.2<-grf(sim,grid=grid,nsim=no.sim,cov.pars=c(1,model$gstat$model[[1]][2,3]),nug=0,mean=0)$data
save(sill.2,file="sill2.Rdata")
rm(sill.2)
gc()





load("sill1.Rdata")

temp1<-sill.1[,1:500]
save(temp1,file="silltemp1.Rdata")
rm(temp1)
gc()

temp2<-sill.1[,501:1000]
save(temp2,file="silltemp2.Rdata")
rm(temp2)
gc()

temp3<-sill.1[,1001:1500]
save(temp3,file="silltemp3.Rdata")
rm(temp3)
gc()

temp4<-sill.1[,1501:2000]
save(temp4,file="silltemp4.Rdata")
rm(temp4)
gc()

temp5<-sill.1[,2001:2500]
save(temp5,file="silltemp5.Rdata")
rm(temp5)
gc()
rm(sill.1)
gc()

load("sill2.Rdata")

temp1<-sill.2[,1:500]
save(temp1,file="sill2temp1.Rdata")
rm(temp1)
gc()

temp2<-sill.2[,501:1000]
save(temp2,file="sill2temp2.Rdata")
rm(temp2)
gc()

temp3<-sill.2[,1001:1500]
save(temp3,file="sill2temp3.Rdata")
rm(temp3)
gc()

temp4<-sill.2[,1501:2000]
save(temp4,file="sill2temp4.Rdata")
rm(temp4)
gc()

temp5<-sill.2[,2001:2500]
save(temp5,file="sill2temp5.Rdata")
rm(temp5)
gc()
rm(sill.2)
gc()



nug.1<-matrix(rnorm(sim*no.sim),ncol=no.sim,nrow=sim)
temp1<-nug.1[,1:500]
save(temp1,file="nugtemp1.Rdata")
rm(temp1)
gc()

temp2<-nug.1[,501:1000]
save(temp2,file="nugtemp2.Rdata")
rm(temp2)
gc()

temp3<-nug.1[,1001:1500]
save(temp3,file="nugtemp3.Rdata")
rm(temp3)
gc()

temp4<-nug.1[,1501:2000]
save(temp4,file="nugtemp4.Rdata")
rm(temp4)
gc()

temp5<-nug.1[,2001:2500]
save(temp5,file="nugtemp5.Rdata")
rm(temp5)
gc()

rm(nug.1)
gc()

nug.2<-matrix(rnorm(sim*no.sim),ncol=no.sim,nrow=sim)
temp1<-nug.2[,1:500]
save(temp1,file="nug2temp1.Rdata")
rm(temp1)
gc()

temp2<-nug.2[,501:1000]
save(temp2,file="nug2temp2.Rdata")
rm(temp2)
gc()

temp3<-nug.2[,1001:1500]
save(temp3,file="nug2temp3.Rdata")
rm(temp3)
gc()

temp4<-nug.2[,1501:2000]
save(temp4,file="nug2temp4.Rdata")
rm(temp4)
gc()

temp5<-nug.2[,2001:2500]
save(temp5,file="nug2temp5.Rdata")
rm(temp5)
gc()

rm(nug.2)
gc()




##############




load("silltemp1.Rdata")
sill.1<-temp1
gc()
load("sill2temp1.Rdata")
sill.2<-temp1
gc()
load('nugtemp1.Rdata')
nug.1<-temp1
gc()
load('nug2temp1.Rdata')
nug.2<-temp1
gc()


sim.tp<-a0[1,1]*nug.1+a0[1,2]*nug.2+a1[1,1]*sill.1+a1[1,2]*sill.2+mean(data$TP)
save(sim.tp,file="simtp1.Rdata")
rm(sim.tp)
gc()
sim.flow<-a0[2,1]*nug.1+a0[2,2]*nug.2+a1[2,1]*sill.1+a1[2,2]*sill.2+mean(data$FLOW)
save(sim.flow,file="simflow1.Rdata")
rm(sim.flow)
gc()
rm(tp.sill)
rm(flow.sill)
rm(nug.tp)
rm(nug.flow)
gc()




load("silltemp2.Rdata")
sill.1<-temp2
gc()
load("sill2temp2.Rdata")
sill.2<-temp2
gc()
load('nugtemp2.Rdata')
nug.1<-temp2
gc()
load('nug2temp2.Rdata')
nug.2<-temp2
gc()


sim.tp<-a0[1,1]*nug.1+a0[1,2]*nug.2+a1[1,1]*sill.1+a1[1,2]*sill.2+mean(data$TP)
save(sim.tp,file="simtp2.Rdata")
rm(sim.tp)
gc()
sim.flow<-a0[2,1]*nug.1+a0[2,2]*nug.2+a1[2,1]*sill.1+a1[2,2]*sill.2+mean(data$FLOW)
save(sim.flow,file="simflow2.Rdata")
rm(sim.flow)
gc()
rm(tp.sill)
rm(flow.sill)
rm(nug.tp)
rm(nug.flow)
gc()



load("silltemp3.Rdata")
sill.1<-temp3
gc()
load("sill2temp3.Rdata")
sill.2<-temp3
gc()
load('nugtemp3.Rdata')
nug.1<-temp3
gc()
load('nug2temp3.Rdata')
nug.2<-temp3
gc()


sim.tp<-a0[1,1]*nug.1+a0[1,2]*nug.2+a1[1,1]*sill.1+a1[1,2]*sill.2+mean(data$TP)
save(sim.tp,file="simtp3.Rdata")
rm(sim.tp)
gc()
sim.flow<-a0[2,1]*nug.1+a0[2,2]*nug.2+a1[2,1]*sill.1+a1[2,2]*sill.2+mean(data$FLOW)
save(sim.flow,file="simflow3.Rdata")
rm(sim.flow)
gc()
rm(tp.sill)
rm(flow.sill)
rm(nug.tp)
rm(nug.flow)
gc()


load("silltemp4.Rdata")
sill.1<-temp4
gc()
load("sill2temp4.Rdata")
sill.2<-temp4
gc()
load('nugtemp4.Rdata')
nug.1<-temp4
gc()
load('nug2temp4.Rdata')
nug.2<-temp4
gc()


sim.tp<-a0[1,1]*nug.1+a0[1,2]*nug.2+a1[1,1]*sill.1+a1[1,2]*sill.2+mean(data$TP)
save(sim.tp,file="simtp4.Rdata")
rm(sim.tp)
gc()
sim.flow<-a0[2,1]*nug.1+a0[2,2]*nug.2+a1[2,1]*sill.1+a1[2,2]*sill.2+mean(data$FLOW)
save(sim.flow,file="simflow4.Rdata")
rm(sim.flow)
gc()
rm(tp.sill)
rm(flow.sill)
rm(nug.tp)
rm(nug.flow)
gc()


load("silltemp5.Rdata")
sill.1<-temp5
gc()
load("sill2temp5.Rdata")
sill.2<-temp5
gc()
load('nugtemp5.Rdata')
nug.1<-temp5
gc()
load('nug2temp5.Rdata')
nug.2<-temp5
gc()


sim.tp<-a0[1,1]*nug.1+a0[1,2]*nug.2+a1[1,1]*sill.1+a1[1,2]*sill.2+mean(data$TP)
save(sim.tp,file="simtp5.Rdata")
rm(sim.tp)
gc()
sim.flow<-a0[2,1]*nug.1+a0[2,2]*nug.2+a1[2,1]*sill.1+a1[2,2]*sill.2+mean(data$FLOW)
save(sim.flow,file="simflow5.Rdata")
rm(sim.flow)
gc()
rm(tp.sill)
rm(flow.sill)
rm(nug.tp)
rm(nug.flow)
gc()













###########



load("simtp1.Rdata")
temp<-sim.tp
rm(sim.tp)
gc()
load("simtp2.Rdata")
temp<-cbind(temp,sim.tp)
rm(sim.tp)
gc()
load("simtp3.Rdata")
temp<-cbind(temp,sim.tp)
rm(sim.tp)
gc()
load("simtp3.Rdata")
temp<-cbind(temp,sim.tp)
rm(sim.tp)
gc()
load("simtp4.Rdata")
sim.tp<-cbind(temp,sim.tp)
rm(temp)
gc()
save(sim.tp,file="simtp.Rdata")
rm(sim.tp)
gc()


load("simflow1.Rdata")
temp<-sim.flow
rm(sim.flow)
gc()
load("simflow1.Rdata")
temp<-cbind(temp,sim.flow)
rm(sim.flow)
gc()
load("simflow3.Rdata")
temp<-cbind(temp,sim.flow)
rm(sim.flow)
gc()
load("simflow4.Rdata")
temp<-cbind(temp,sim.flow)
rm(sim.flow)
gc()
load("simflow5.Rdata")
sim.flow<-cbind(temp,sim.flow)
rm(temp)
gc()
save(sim.flow,file="simflow.Rdata")
gc()

########################################################################
########BACKTRANSFORM THE SIMULATED DATA USING GEOR PACKAGE.############
########################################################################
###Calculate the variance
load("simflow.Rdata")
sim.flow.var<-apply(sim.flow,1,var)
###Create a list for the results
simulatedFlow<-matrix(NA,ncol=2500,nrow=length(sim.flow.var))
for(i in 1:2500){
simulatedFlow[,i]<-backtransform.moments(lambda=flow.lambda,mean=sim.flow[,i],sim.flow.var)$mean
print(i)
}
rm(sim.flow)
rm(sim.flow.var)
gc()
save(simulatedFlow,file="simulatedFlow.Rdata")
rm(simulatedFlow)
gc()
###Calculate the variance
load("simtp.Rdata")
sim.tp.var<-apply(sim.tp,1,var)
###Create a list for the results
simulatedTP<-matrix(NA,ncol=2500,nrow=length(sim.tp.var))
for(i in 1:2500){
simulatedTP[,i]<-backtransform.moments(lambda=tp.lambda,mean=sim.tp[,i],sim.tp.var)$mean
}
rm(sim.tp)
rm(sim.tp.var)
gc()
save(simulatedTP,file="simulatedTP.Rdata")
rm(simulatedTP)
gc()
##############################################################################
#### Simulate sampling of the data using event-based and routine sampling ####
##############################################################################


##########################
#### Routine sampling ####
##########################


#### Create a seq of dates for one year with each day repeated 24 times
date<-rep(seq.Date(as.Date("1991-01-01"),as.Date("1991-12-31"),by="days"),each=24)
#Extract the numeric day of each date
date<-as.numeric(format(date,"%d"))
#Join the day number with the hour of the day   
date<-paste(date,rep(seq(0:23),365)-1,sep=".")
#Then repeat the sequence above (hourly seq of one year) for 20 years
date<-rep(date,20)
###Create a custom time to represtent when the routine sample was taken.
routineTime<-(1:length(date))/24
routineTime<-routineTime[date==15.12]
###Subset the simulated data on the 15th of each month at miday.
load("simulatedTP.Rdata")
routineTP<-simulatedTP[date==15.12,]
save(routineTP,file="routineTP.Rdata")
rm(simulatedTP)
gc()
load("simulatedFlow.Rdata")
routineFlow<-simulatedFlow[date==15.12,]
save(routineFlow,file="routineFlow.Rdata")
rm(simulatedFlow)
gc()

#############################################
##--------- Event based sampling ----------##
#############################################
load("simulatedFlow.Rdata")
load("simulatedTP.Rdata")
eventData<-list()
###Need to convert the Discharge (ML/day) to (m/day), using a rating curve
ratingCurve<-read.csv("rating_curve.csv",header=T)
###Fit a model to the rating curve.
ratingCurveModel<-smooth.spline(ratingCurve[,2],ratingCurve[,1])
###Create a vector of the amount of hours between each sample
sampleHour<-c(0,3,3,3,6,6,6,6,12,12,12,12)

for(realisation in 1:2500){
	###Create a height object from the above model
	streamHeight<-predict(ratingCurveModel,simulatedFlow[,realisation]*(24))$y
	###Create a vector of 
	###Create Variables for time of sample, discharge and TP.
	tempTime<-NA
	tempFlow<-NA
	tempTP<-NA
	###Set the while loop iterator to 2
	i=1
	###Now sample the discharge if there is an event.
	while(i<length(streamHeight)-81){
		i=i+1
		##Check to see if the previous hour of stream height was less 1 m and if the current stream height is more than 1 m.
		if(streamHeight[i-1]<1&&streamHeight[i]>1){
			##Check to see that the rate of rise is enough to classify this as an event.
			if(streamHeight[i+1]-streamHeight[i]>0.04*2){


	##----Now we are within an event.----##


				##Set the iterator of the event to 1. This iterator will cycle through the sampleHour
				j=1
				while(streamHeight[i]>=1&&j<12){#stream height must be more than 1 m high, and cannot take more than 12 samples.
					i=i+sampleHour[j] #This makes sure i will be at the right hour of time when we leave the event.
					####Now save the associated values at each time of the event.
					tempFlow<-rbind(tempFlow,simulatedFlow[i,realisation])
					tempTP<-rbind(tempTP,simulatedTP[i,realisation])
					tempTime<-rbind(tempTime,i)
					j=j+1
				}
			}
		}	
	}

	###Combine the vectors to one data.frame and remove the first obs, as they are NA's.
	eventSampled<-data.frame(Time=(tempTime[-1]/24),Flow=tempFlow[-1],TP=tempTP[-1])
	###From here i am trying to join the routine and the event	
	tempfinal<-data.frame(Time=routineTime,Flow=routineFlow[,realisation],TP=routineTP[,realisation])
	tempfinal<-rbind(eventSampled,tempfinal)
	tempfinal<-tempfinal[order(tempfinal$Time),]
	eventData[[realisation]]<-tempfinal
	print(realisation)
}
save(eventData,file="eventData.Rdata")

###Lets see how many samples where collected for each realisation.
temp<-length(eventData[[1]][,1])
for(i in 2:1410){
temp<-rbind(temp,length(eventData[[i]][,1]))
if(temp[i]==2793)print(i)
}
hist(temp)
####So now there is both routine and event-based sampling completed - all that has to be done is to make likfit models of the whole lot.

################################
###Create the likfit objects
setwd("~/Documents/code/Simulations/")
load("eventData.Rdata")
library(geoR)
eventlikfit<-list()

ini.=matrix(NA,ncol=2,nrow=13)
ini.[,1]=c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)
ini.[,2]=c(2,3,4,5,6,7,8,10,20,30,50,70,100)
j=1
for(i in 2251:2500){


	lambda<-boxcox.fit(eventData[[i]]$TP)$lambda
	flowlambda<-boxcox.fit(eventData[[i]]$Flow)$lambda
#eventData$flow<-BCtransform(eventData$flow,flowlambda)$data
##Create a geodata object of the data.
#BCtransform(eventData[[i]]$Flow,flowlambda)$data
	eventOKData<-data.frame(x=eventData[[i]]$Time,y=1,tp=eventData[[i]]$TP,flow=BCtransform(eventData[[i]]$Flow,flowlambda)$data)
	eventOKData<-eventOKData[!duplicated(eventOKData),]
	eventGeo<-as.geodata(eventOKData,coords.col=1:2,data.col=3,covar.col=4)
###Plot a variogram of the data, using lambda to transform the data
#plot(variog(eventGeo,lambda=lambda,max.dist=30,trend=~flow))
###Use likfit to optimise the variogram model.
		eventlikfit[[j]]<-likfit(eventGeo,ini=ini.,fix.nugget = F,fix.lambda=T, lik.method = "REML",lambda=lambda,trend=trend.spatial(~flow))
#lines(eventlikfit)
j=j+1
print(j)
}
save(eventlikfit,file="eventlikfit_part10.Rdata")
rm(eventlikfit)
gc()

load("simulatedFlow.Rdata")
temp<-simulatedFlow[,1:250]
save(temp,file="simflow1.Rdata")
rm(temp)
gc()
temp<-simulatedFlow[,251:500]
save(temp,file="simflow2.Rdata")
rm(temp)
gc()
temp<-simulatedFlow[,501:750]
save(temp,file="simflow3.Rdata")
rm(temp)
gc()
temp<-simulatedFlow[,751:1000]
save(temp,file="simflow4.Rdata")
rm(temp)
gc()
temp<-simulatedFlow[,1001:1250]
save(temp,file="simflow5.Rdata")
rm(temp)
gc()
temp<-simulatedFlow[,1251:1500]
save(temp,file="simflow6.Rdata")
rm(temp)
gc()
temp<-simulatedFlow[,1501:1750]
save(temp,file="simflow7.Rdata")
rm(temp)
gc()
temp<-simulatedFlow[,1751:2000]
save(temp,file="simflow8.Rdata")
rm(temp)
gc()
temp<-simulatedFlow[,2001:2250]
save(temp,file="simflow9.Rdata")
rm(temp)
gc()
temp<-simulatedFlow[,2251:2500]
save(temp,file="simflow10.Rdata")
rm(temp)
gc()



load("simulatedTP.Rdata")
temp<-simulatedTP[,1:250]
save(temp,file="simtp1.Rdata")
rm(temp)
gc()
temp<-simulatedTP[,251:500]
save(temp,file="simtp2.Rdata")
rm(temp)
gc()
temp<-simulatedTP[,501:750]
save(temp,file="simtp3.Rdata")
rm(temp)
gc()
temp<-simulatedTP[,751:1000]
save(temp,file="simtp4.Rdata")
rm(temp)
gc()
temp<-simulatedTP[,1001:1250]
save(temp,file="simtp5.Rdata")
rm(temp)
gc()
temp<-simulatedTP[,1251:1500]
save(temp,file="simtp6.Rdata")
rm(temp)
gc()
temp<-simulatedTP[,1501:1750]
save(temp,file="simtp7.Rdata")
rm(temp)
gc()
temp<-simulatedTP[,1751:2000]
save(temp,file="simtp8.Rdata")
rm(temp)
gc()
temp<-simulatedTP[,2001:2250]
save(temp,file="simtp9.Rdata")
rm(temp)
gc()
temp<-simulatedTP[,2251:2500]
save(temp,file="simtp10.Rdata")
rm(temp)
gc()




####Get the distances
setwd("~/Documents/code/Simulations/")
library(geoR)
date<-rep(seq.Date(as.Date("1991-01-01"),as.Date("1991-12-31"),by="days"),each=24)
#Extract the numeric day of each date
date<-as.numeric(format(date,"%d"))
#Join the day number with the hour of the day   
date<-paste(date,rep(seq(0:23),365)-1,sep=".")
#Then repeat the sequence above (hourly seq of one year) for 20 years
date<-rep(date,20)
###Create a custom time to represtent when the routine sample was taken.
routineTime<-(1:length(date))/24
routineTime<-routineTime[date==15.12]

##############
load("routineFlow.Rdata")
load("routineTP.Rdata")

ini.=matrix(NA,ncol=2,nrow=6)
ini.[,1]=c(0.01,0.01,0.01,0.01,0.01,0.01)
ini.[,2]=c(100,20,30,10,150,300)
routinelikfit<-list()
for(i in 1:2500){
	lambda<-boxcox.fit(routineTP[,i])$lambda
	flowlambda<-boxcox.fit(routineFlow[,i])$lambda

	routineOKData<-data.frame(x=routineTime,y=1,tp=routineTP[,i],flow=BCtransform(routineFlow[,i],flowlambda)$data)
	routineGeo<-as.geodata(routineOKData,coords.col=1:2,data.col=3,covar.col=4)

	routinelikfit[[i]]<-likfit(routineGeo,ini=ini.,fix.nugget = F,fix.lambda=T, lik.method = "REML",lambda=lambda,trend=trend.spatial(~flow))
	print(i)
}

save(routinelikfit,file="routinelikfit.Rdata")
###Extract the AIC values for the non spatial and the spatial models
spatial<-NA
ns.spatial<-NA
range<-NA
real<-NA
for(i in 1:2500){
	spatial<-rbind(spatial,summary(routinelikfit[[i]])$likelihood$AIC)
	ns.spatial<-rbind(ns.spatial,summary(routinelikfit[[i]])$ns.likelihood$AIC)
	if(i>1){
		if(spatial[i]<ns.spatial[i]){
			real<-rbind(real,i)
			range<-rbind(range,summary(routinelikfit[[i]])$practicalRange/3)
		}
	}
}

real<-as.numeric(real[-1,])
range<-as.numeric(range[-1,])
ns.spatial<-as.numeric(ns.spatial[-1,])
spatial<-as.numeric(spatial[-1,])
diff<-spatial-ns.spatial
par(mfrow=c(2,2))
hist(ns.spatial)
hist(spatial)
hist(diff)
hist(range)

routinelikfitproperties<-list()
routinelikfitproperties[["spatialAIC"]]<-spatial
routinelikfitproperties[["non.spatialAIC"]]<-ns.spatial
routinelikfitproperties[["AIC.difference"]]<-diff
routinelikfitproperties[["spatial.Realisations"]]<-real
routinelikfitproperties[["spatial.range"]]<-range
save(routinelikfitproperties,file="routinelikfitproperties.Rdata")
#

#Find out how many realisations had a spatial model- with a range less than 365
length(real)


######################################################
##### now look at how close the simulations are ######
##### by doing lmcr's of the simulated data ##########
######################################################

seq<-1:2500
n<-100
random.cols<-sample(x=seq,size=n,replace=FALSE)

setwd("~/Documents/code/Simulations/")
load("simulatedFlow.Rdata")
smallFlow<-simulatedFlow[,random.cols]
rm(simulatedFlow)
gc()
load("simulatedTP.Rdata")
smallTP<-simulatedTP[,random.cols]
rm(simulatedTP)
gc()



### Load the gstat library
library(gstat);library(TeachingDemos);library(geoR)
sim=24*365*20
date=(seq(1,sim)-1)*0.04166667
small.gstat<-list()
for(i in 1:100){
	data <- data.frame(X=date,Y=1,TP=smallTP[,i],FLOW=smallFlow[,i])
	flow.lambda <- boxcox.fit(data$FLOW+0.001)$lambda
	tp.lambda <- boxcox.fit(data$TP+0.001)$lambda
### Transform the data the boxcox formula.
	data$TP <- ((data$TP+0.001)^tp.lambda-1)/tp.lambda
	data$FLOW <- ((data$FLOW+0.001)^flow.lambda-1)/flow.lambda
### Create gstat object of the data
	spdf <- SpatialPointsDataFrame(data[,1:2],data)
### Create A gstat object with TP and flow.
	g = gstat(NULL, "TP", TP ~ 1,spdf,maxdist=200,dummy=TRUE,beta=mean(data$TP))
	g = gstat(g, "FLOW", FLOW ~ 1,spdf,maxdist=200,dummy=TRUE,beta=mean(data$FLOW))


### Create cross and auto variograms of the gstat object, using width to set the bins.
	small.gstat[[i]] = variogram(g,cutoff=150,width=10)
print(i)
}
save(small.gstat,file="sim_subset_gstat_validation.Rdata")
###Now to have a look at the results
results<-matrix(rep(NA,45*100),ncol=100)
for(i in 1:100){
results[,i]<-small.gstat[[i]]$gamma
}
mean<-apply(results,1,summary)[4,]

demo<-small.gstat[[1]]
demo$gamma<-mean
plot(demo,model$gstat)



