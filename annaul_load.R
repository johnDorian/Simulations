###This script is designed to get the annual load, concentration and flow from the actual and predicted data sets. It is also designed to determine how many hours the threshold is exceeded (0.2). 

###TODO: Re-create the SCA method of giving a percentage of samples (dry and wet and seperated) above the threshold. Then compare this with actual results, and the two predicted results.


realisation=0 #An iterator for 1:2500 realisations
###Now some matrixs for the results - 2500 realisations, and 20 years in each realisation
eventTP<-matrix(nrow=20,ncol=2500)
eventmodTP<-matrix(nrow=20,ncol=2500)
routine.TP<-matrix(nrow=20,ncol=2500)
routinemodTP<-matrix(nrow=20,ncol=2500)
actualTP<-matrix(nrow=20,ncol=2500)
actualFlow<-matrix(nrow=20,ncol=2500)
#Firstly we will get the annual load for the predicted things
for(part in 1:10){
###Load the simulated TP and Flow for the next 250 realisations.
load(paste('~/Documents/code/Simulations/data/backtransformed_simulations/parts/tp/simulatedTP',part,'.Rdata',sep=""))
load(paste('~/Documents/code/Simulations/data/backtransformed_simulations/parts/flow/simulatedFlow',part,'.Rdata',sep=""))
###Load the sampled data.
load(paste('~/Documents/code/Simulations/data/simulated_sampled_data/event/eventDatapart',part,'.Rdata',sep=""))
load(paste('~/Documents/code/Simulations/data/simulated_sampled_data/routine/routineTPpart',part,'.Rdata',sep=""))
	for(subpart in 1:250){
		##Get the actual results
		realisation<-realisation+1 # just a counter from 1 to 2500 - for each realisation
		###Seperate the realisation into a mtrix of years and get the sum of each year.
		actualTP[,realisation]<-colSums(matrix(simulatedTP[,subpart],ncol=20)) 
		actualFlow[,realisation]<-colSums(matrix(simulatedFlow[,subpart],ncol=20))
		##Load the predicted values for this realisation
		load(paste('~/Documents/code/Simulations/predicted/event/part',part,'/krigtp_part',part,'_subpart',subpart,'.Rdata',sep=""))
		eventkrig<-krig.tp
		load(paste('~/Documents/code/Simulations/predicted/routine/part',part,'/krigtp_',part,'_subpart',subpart,'.Rdata',sep=""))
		###Get the max tp value
		max.E.TP<-max(eventData[[subpart]]$TP)*2
		max.R.TP<-max(routineTP[,subpart])*2
		##Go through this relisation and get the sum of each year.
		for(year in 1:20){
			###Calculate the annual load of each sampling method. Each sampling method has two load estimations, as there is now a limit for TP maximum which is 2*maximum sampled TP,  using the assoicated method.
			#First method is event kriging
			temp<-eventkrig
			temp[[year]][[1]][eventkrig[[year]][[1]]>max.E.TP]=max.E.TP
			eventmodTP[year,realisation]<-sum(temp[[year]][[1]])
			eventTP[year,realisation]<-sum(eventkrig[[year]][[1]])
			#Second method is routine
			temp<-krig.tp
			temp[[year]][[1]][krig.tp[[year]][[1]]>max.R.TP]=max.R.TP
			routinemodTP[year,realisation] <- sum(temp[[year]][[1]])
			routine.TP[year,realisation]<-sum(krig.tp[[year]][[1]])
		}
	print(realisation)
	}

}
##Now calculate the actual loads...
actualLoad<-actualTP*actualFlow
eventLoad<-eventTP*actualFlow
eventmodLoad<-eventmodTP*actualFlow
routineLoad<-routine.TP*actualFlow
routinemodLoad<-routinemodTP*actualFlow
#Now save all this so we don't have to do this all over again.
all.years <- list(ActLoad=actualLoad,EvLoad=eventLoad,EvModLoad=eventmodLoad,RoLoad=routineLoad,RoModLoad=routinemodLoad,ActFlow=actualFlow,ActTP=actualTP,roTP=routine.TP,roModTP=routinemodTP,evTP=eventTP,evModTP=eventmodTP)
save(all.years,file="~/Documents/code/Simulations/AnnualLoads.Rdata")



###See what the new work does to the results.
#Set the color system - for transparent colors.
gr <- rgb(green=255,red=0,blue=0,alpha=100,max=255)
rd <- rgb(green=0,red=255,blue=0,alpha=100,max=255)
bl <- rgb(green=0,red=0,blue=0,alpha=100,max=255)
bu <- rgb(green=0,red=0,blue=255,alpha=100,max=255)
#Close any open graphic windows.
dev.off()
###Problem realisation - the one with crazy TP values.
realisation=1871
###Random realisation
#realisation=sample(2500,1)

##Plot the four different annual loads using the color system from above.
jpeg("Problem_realisation_fixed.jpg")
plot(eventLoad[,realisation],col=rd,pch=16,cex=0.8,xlab="Year of realisation",ylab="Annual load (mg)")
points(actualLoad[,realisation],pch=16,col=bl,cex=0.8)
points(routineLoad[,realisation],col=bu,pch=16,cex=0.8)
points(eventmodLoad[,realisation],col=gr,pch=16,cex=0.8)
#Add a nice pretty legend for others to know what is going on.
legend("topright",c("actual","event","event modified",'routine'),col=c(bl,rd,gr,bu),pch=16)
dev.off()

###Now look at how much this modification method has changed the results
jpeg("How_similar.jpg")
##Get the min and maximum values of both the x and y variables
ylim=c(min(c(eventLoad,eventmodLoad)),max(c(eventLoad,eventmodLoad)))
#Set the above limits for the x axis limits
xlim=ylim
#Plot the two annual loads using the limits and the transparent colours, and the axis limits
plot(eventmodLoad,eventLoad,ylim=ylim,xlim=xlim,col=bl,pch=16,cex=0.6)
#Now add a 45 degree line to the graph.
abline(0,1,lty=2,col="red")
##How many really are different - using as.numeric for binary values 1 and 0, true and false, therefore you can use sum to find the number of similar.
legend("bottomright",paste("Years changed:",sum(as.numeric(eventmodLoad!=eventLoad))))
dev.off()


ylim=c(min(c(routineLoad,routinemodLoad)),max(c(routineLoad,routinemodLoad)))
#Set the above limits for the x axis limits
xlim=ylim
#Plot the two annual loads using the limits and the transparent colours, and the axis limits
plot(routinemodLoad,routineLoad,ylim=ylim,xlim=xlim,col=bl,pch=16,cex=0.6)
#Now add a 45 degree line to the graph.
abline(0,1,lty=2,col="red")





###############################
#####This next section will look at how many hours in the year the tp is over the threshold of 0.2
###############################
for(part in 1:10){
	load(paste("~/Documents/code/Simulations/data/backtransformed_simulations/parts/tp/simulatedTP",part,".Rdata",sep=""))
	for(subpart in 1:250){
		real.tp<-matrix(simulatedTP[,subpart],ncol=20)



###The results go in here
hours.above<-list()
###The rmse results go in here
rmse.results<-matrix(NA,ncol=2,nrow=2500)

#The rmse function
rmse <- function(obs, pred) sqrt(mean((obs-pred)^2))
#The iterator for the each realisation
iter <- 0

for(part in 1:10){

	###Load the sampled data.
	load(paste('~/Documents/code/Simulations/data/simulated_sampled_data/event/eventDatapart',part,'.Rdata',sep=""))
	load(paste('~/Documents/code/Simulations/data/simulated_sampled_data/routine/routineTPpart',part,'.Rdata',sep=""))
	##Load the TP data
	load(paste("data/backtransformed_simulations/parts/tp/simulatedTP",part,".Rdata",sep=""))

	##Cycle through each realisation of each part
	for(subpart in 1:250){
		##Increase the iterator
		iter <- iter + 1
		##Put this relisation of TP into a nice matrix for each year of the realisation
		real.tp<-matrix(simulatedTP[,subpart],ncol=20)
		###Load the realisation of event prediction
		load(paste("predicted/event/part",part,"/krigtp_part",part,"_subpart",subpart,".Rdata",sep=""))
		###Create a 
		event.tp<-matrix(NA,ncol=20,nrow=8760)
		###Get the maximum value of the sampled data of the two sampling methods
		max.E.TP<-max(eventData[[subpart]]$TP)*2
		max.R.TP<-max(routineTP[,subpart])*2		
		for(i in 1:20){
			#Change any predicted values to that of the maximum limit, based on the 2*maximum of the sampled data
			temp<-krig.tp
			temp[[i]][[1]][krig.tp[[i]][[1]]>max.E.TP]=max.E.TP
			event.tp[,i]<-temp[[i]][[1]]
		}
		routine.tp<-matrix(NA,ncol=20,nrow=8760)
		load(paste("predicted/routine/part",part,"/krigtp_",part,"_subpart",subpart,".Rdata",sep=""))		
		for(i in 1:20){
			#Change any predicted values to that of the maximum limit, based on the 2*maximum of the sampled data
			temp<-krig.tp
			temp[[i]][[1]][krig.tp[[i]][[1]]>max.R.TP]=max.R.TP
			routine.tp[,i]<-temp[[i]][[1]]
		}
		days=data.frame(actual=rep(NA,20),event=NA,routine=NA)
		for(i in 1:20){
			days[i,1]<-length(real.tp[real.tp[,i]>0.2,1])
			days[i,2]<-length(event.tp[event.tp[,i]>0.2,1])
			days[i,3]<-length(routine.tp[routine.tp[,i]>0.2,1])
		}
		hours.above[[(((part-1)*250)+subpart)]]<-days
		rmse.results[(((part-1)*250)+subpart),1:2]<-c(rmse(days[,1],days[,2]),rmse(days[,1],days[,3]))
	}
	print(part)
}
save(hours.above,file="~/Documents/code/Simulations/hoursAbove.Rdata")
save(rmse.results,file="~/Documents/code/Simulations/hoursAboveRMSE.Rdata")

#######################################################
####Section 3: Recreation of the SCA exceedance tables
#######################################################

###Results matrix
eTotal<-matrix(nrow=20,ncol=2500)
eOver<-matrix(nrow=20,ncol=2500)
rOver<-matrix(nrow=20,ncol=2500)

##Get the start and ends of each month
time<-matrix(seq(1/24,(365*20),1/24),ncol=20)
time<-time[c(1,nrow(time)),]

realisation=0
for(part in 1:10){

	###Load the sampled data fore both sampling methods
	load(paste('~/Documents/code/Simulations/data/simulated_sampled_data/event/eventDatapart',part,'.Rdata',sep=""))
	load(paste('~/Documents/code/Simulations/data/simulated_sampled_data/routine/routineTPpart',part,'.Rdata',sep=""))


	for(subpart in 1:250){
		realisation=realisation+1
#########First section 
		rOver[,realisation]<-colSums(matrix(as.numeric(routineTP[,subpart]>0.2),ncol=20))


###Now for the other method.
###Subset the year of observations.
		for(year in 1:20){
			temp<-eventData[[subpart]][eventData[[subpart]][[1]]>=time[1,year]&eventData[[subpart]][[1]]<time[2,year],]
			eOver[year,realisation]<-sum(as.numeric(temp$TP>0.2))
			eTotal[year,realisation]<-length(temp$TP)
		}
	}
}
###Save the results
thresholdSummary<-list(eventOver=eOver,eventTotal=eTotal,routineOver=rOver)
save(thresholdSummary,file="~/Documents/code/Simulations/thresholdSummary.Rdata")

hoursAbove<-do.call(cbind,hours.above)
hoursActual<-as.numeric(do.call(rbind,hoursAbove[,seq(1,7500,3)]))
hoursEvent<-as.numeric(do.call(rbind,hoursAbove[,seq(2,7500,3)]))
hoursRoutine<-as.numeric(do.call(rbind,hoursAbove[,seq(3,7500,3)]))

plot(hoursActual/(365*24),hoursEvent/(365*24))
plot(hoursActual/(365*24),hoursRoutine/(365*24))






