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


load("~/Documents/code/Simulations/data/lmcr_model/lmcr.Rdata")


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
source('~/Documents/code/Simulations/functions/fiting.R')
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


##############REDO EVERYTHING FROM HERE 9/08/2010...


load("~/Documents/code/Simulations/data/simulated_raw/sills/sill1.Rdata")
start=seq(1,2251,250)
end=seq(250,2500,250)

for(i in 1:length(start)){
	temp<-sill.1[,start[i]:end[i]]
	save(temp,file=paste("~/Documents/code/Simulations/data/simulated_raw/sills/sill1temp",i,".Rdata",sep=""))
	rm(temp)
	gc()
}
rm(sill.1)
gc()

load("~/Documents/code/Simulations/data/simulated_raw/sills/sill2.Rdata")
for(i in 1:length(start)){
	temp<-sill.2[,start[i]:end[i]]
	save(temp,file=paste("~/Documents/code/Simulations/data/simulated_raw/sills/sill2temp",i,".Rdata",sep=""))
	rm(temp)
	gc()
}
rm(sill.2)
gc()

sim=24*365*20
no.sim=2500
nug.1<-rnorm(sim*no.sim)

start=seq(1,((length(nug.1)-sim)+1),sim*250)
end=seq(sim*250,length(nug.1),sim*250)


for(i in 1:length(start)){
	temp<-matrix(nug.1[start[i]:end[i]],ncol=250)
	save(temp,file=paste("~/Documents/code/Simulations/data/simulated_raw/nuggets/nug1temp",i,".Rdata",sep=""))
	rm(temp)
	gc()
}
rm(nug.1)
gc()

nug.2<-rnorm(sim*no.sim)
for(i in 1:length(start)){
	temp<-matrix(nug.2[start[i]:end[i]],ncol=250)
	save(temp,file=paste("~/Documents/code/Simulations/data/simulated_raw/nuggets/nug2temp",i,".Rdata",sep=""))
	rm(temp)
	gc()
}
rm(nug.2)
gc()

############# Create the simulations...
for(i in 1:10){
	load(paste("~/Documents/code/Simulations/data/simulated_raw/sills/sill1temp",i,".Rdata",sep=""))
	sill.1<-temp
	gc()
	load(paste("~/Documents/code/Simulations/data/simulated_raw/sills/sill2temp",i,".Rdata",sep=""))
	sill.2<-temp
	gc()
	load(paste("~/Documents/code/Simulations/data/simulated_raw/nuggets/nug1temp",i,".Rdata",sep=""))
	nug.1<-temp
	gc()
	load(paste("~/Documents/code/Simulations/data/simulated_raw/nuggets/nug2temp",i,".Rdata",sep=""))
	nug.2<-temp
	gc()

	sim.tp<-a0[1,1]*nug.1+a0[1,2]*nug.2+a1[1,1]*sill.1+a1[1,2]*sill.2+mean(data$TP)
	save(sim.tp,file=paste("~/Documents/code/Simulations/data/simulated_raw/simtp/simtp",i,".Rdata",sep=""))
	rm(sim.tp)
	gc()
	sim.flow<-a0[2,1]*nug.1+a0[2,2]*nug.2+a1[2,1]*sill.1+a1[2,2]*sill.2+mean(data$FLOW)
	save(sim.flow,file=paste("~/Documents/code/Simulations/data/simulated_raw/simflow/simflow",i,".Rdata",sep=""))
	rm(sim.flow)

	rm(sill.1)
	rm(sill.2)
	rm(nug.1)
	rm(nug.2)
	gc()
}

########################################################################
########BACKTRANSFORM THE SIMULATED DATA USING GEOR PACKAGE.############
########################################################################

###First section loads all the simulated flow and gets the variance of it.
load(paste("~/Documents/code/Simulations/data/simulated_raw/simflow/simflow1.Rdata",sep=""))
all.sim.flow1<-sim.flow[1:87600,]
for(i in 2:10){
	load(paste("~/Documents/code/Simulations/data/simulated_raw/simflow/simflow",i,".Rdata",sep=""))
	all.sim.flow1<-cbind(all.sim.flow1,sim.flow[1:87600,])
}
rm(sim.flow)
gc()
sim.flow.var1<-apply(all.sim.flow1,1,var)
rm(all.sim.flow1)
gc()

load(paste("~/Documents/code/Simulations/data/simulated_raw/simflow/simflow1.Rdata",sep=""))
all.sim.flow2<-sim.flow[87601:175200,]
for(i in 2:10){
	load(paste("~/Documents/code/Simulations/data/simulated_raw/simflow/simflow",i,".Rdata",sep=""))
	all.sim.flow2<-cbind(all.sim.flow2,sim.flow[87601:175200,])
}
rm(sim.flow)
gc()
sim.flow.var2<-apply(all.sim.flow2,1,var)
rm(all.sim.flow2)
gc()
temp<-cbind(sim.flow.var1,sim.flow.var2)
sim.flow.var<-matrix(temp,ncol=1)
rm(temp)
rm(sim.flow.var1)
rm(sim.flow.var2)
gc()

library(geoR)



####Now with the variance we can back transform all the flow
for(i in 1:10){
	simulatedFlow<-matrix(NA,ncol=250,nrow=(365*24*20))
	load(paste("~/Documents/code/Simulations/data/simulated_raw/simflow/simflow",i,".Rdata",sep=""))
	for(real. in 1:250){
		simulatedFlow[,real.]<-backtransform.moments(lambda=flow.lambda,mean=sim.flow[,real.],sim.flow.var)$mean
	}
	save(simulatedFlow,file=paste("~/Documents/code/Simulations/data/backtransformed_simulations/parts/flow/simulatedFlow",i,".Rdata",sep=""))
	rm(simulatedFlow)
	gc()
}
###Now do it all over again for tp
load(paste("~/Documents/code/Simulations/data/simulated_raw/simtp/simtp1.Rdata",sep=""))
all.sim.tp1<-sim.tp[1:87600,]
for(i in 2:10){
	load(paste("~/Documents/code/Simulations/data/simulated_raw/simtp/simtp",i,".Rdata",sep=""))
	all.sim.tp1<-cbind(all.sim.tp1,sim.tp[1:87600,])
}
rm(sim.tp)
gc()
sim.tp.var1<-apply(all.sim.tp1,1,var)
rm(all.sim.tp1)
gc()

load(paste("~/Documents/code/Simulations/data/simulated_raw/simtp/simtp1.Rdata",sep=""))
all.sim.tp2<-sim.tp[87601:175200,]
for(i in 2:10){
	load(paste("~/Documents/code/Simulations/data/simulated_raw/simtp/simtp",i,".Rdata",sep=""))
	all.sim.tp2<-cbind(all.sim.tp2,sim.tp[87601:175200,])
}
rm(sim.tp)
gc()
sim.tp.var2<-apply(all.sim.tp2,1,var)
rm(all.sim.tp2)
gc()
temp<-cbind(sim.tp.var1,sim.tp.var2)
sim.tp.var<-matrix(temp,ncol=1)
rm(temp)
rm(sim.tp.var1)
rm(sim.tp.var2)
gc()


for(i in 1:10){
	simulatedTP<-matrix(NA,ncol=250,nrow=(365*24*20))
	load(paste("~/Documents/code/Simulations/data/simulated_raw/simtp/simtp",i,".Rdata",sep=""))
	for(real. in 1:250){
		simulatedTP[,real.]<-backtransform.moments(lambda=tp.lambda,mean=sim.tp[,real.],sim.tp.var)$mean
	}
	save(simulatedTP,file=paste("~/Documents/code/Simulations/data/backtransformed_simulations/parts/tp/simulatedTP",i,".Rdata",sep=""))
	rm(simulatedTP)
	gc()
}



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

for(i in 1:10){
###Subset the simulated data on the 15th of each month at miday.
	load(paste("~/Documents/code/Simulations/data/backtransformed_simulations/parts/tp/simulatedTP",i,".Rdata",sep=""))
	routineTP<-simulatedTP[date==15.12,]
	save(routineFlow,file=paste("~/Documents/code/Simulations/data/simulated_sampled_data/routineTPpart",i,".Rdata",sep=""))
	rm(simulatedTP)
	gc()
	load(paste("~/Documents/code/Simulations/data/backtransformed_simulations/parts/flow/simulatedFlow",i,".Rdata",sep=""))
	routineFlow<-simulatedFlow[date==15.12,]
	save(routineFlow,file=paste("~/Documents/code/Simulations/data/simulated_sampled_data/routineFlowpart",i,".Rdata",sep=""))
	rm(simulatedFlow)
	gc()
}

#############################################
##--------- Event based sampling ----------##
#############################################

eventData<-list()
###Need to convert the Discharge (ML/day) to (m/day), using a rating curve
ratingCurve<-read.csv("~/Documents/code/Simulations/data/raw_data/rating_curve.csv",header=T)
###Fit a model to the rating curve.
ratingCurveModel<-smooth.spline(ratingCurve[,2],ratingCurve[,1])
###Create a vector of the amount of hours between each sample
sampleHour<-c(0,3,3,3,6,6,6,6,12,12,12,12)

for(file in 1:10){
	load(paste("~/Documents/code/Simulations/data/backtransformed_simulations/parts/flow/simulatedFlow",file,".Rdata",sep=""))
	load(paste("~/Documents/code/Simulations/data/backtransformed_simulations/parts/tp/simulatedTP",file,".Rdata",sep=""))
	load(paste("~/Documents/code/Simulations/data/simulated_sampled_data/routine/routineFlowpart",file,".Rdata",sep=""))
	load(paste("~/Documents/code/Simulations/data/simulated_sampled_data/routine/routineTPpart",file,".Rdata",sep=""))

	for(realisation in 1:250){
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
	save(eventData,file=paste("~/Documents/code/Simulations/data/simulated_sampled_data/event/eventDatapart",file,".Rdata",sep=""))
}

###This section doesn't work as all the samples are seperated now.
###Lets see how many samples where collected for each realisation.
#temp<-length(eventData[[1]][,1])
#for(i in 2:1410){
#temp<-rbind(temp,length(eventData[[i]][,1]))
#if(temp[i]==2793)print(i)
#}
#hist(temp)
####So now there is both routine and event-based sampling completed - all that has to be done is to make likfit models of the whole lot.

###This section is old needs rewritting - 9/Aug/2010.
################################
###Create the likfit objects####
################################
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



setwd("~/Documents/code/Simulations/data/backtransformed_simulations/complete")
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
n<-250
random.cols<-sample(x=seq,size=n,replace=FALSE)
sort(random.cols)
breaks<-seq(0,2500,by=250)
real.<-list()
for(i in 1:(length(breaks)-1)){
	real.[[i]]<-sort(random.cols[random.cols<breaks[i+1]& random.cols>breaks[i]])
}
	




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



