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


eventData<-list()
###Need to convert the Discharge (ML/day) to (m/day), using a rating curve
ratingCurve<-read.csv("~/Documents/code/Simulations/data/raw_data/rating_curve.csv",header=T)
###Fit a model to the rating curve.
ratingCurveModel<-smooth.spline(ratingCurve[,2],ratingCurve[,1])
###Create a vector of the amount of hours between each sample
sampleHour<-c(0,3,3,3,6,6,6,6,12,12,12,12)

for(file in 1:1){
	load(paste("~/Documents/code/Simulations/data/backtransformed_simulations/parts/flow/simulatedFlow",file,".Rdata",sep=""))
	load(paste("~/Documents/code/Simulations/data/backtransformed_simulations/parts/tp/simulatedTP",file,".Rdata",sep=""))
	load(paste("~/Documents/code/Simulations/data/simulated_sampled_data/routine/routineFlowpart",file,".Rdata",sep=""))
	load(paste("~/Documents/code/Simulations/data/simulated_sampled_data/routine/routineTPpart",file,".Rdata",sep=""))

	for(realisation in 1:1){
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
	#save(eventData,file=paste("~/Documents/code/Simulations/data/simulated_sampled_data/event/eventDatapart",file,".Rdata",sep=""))
}


par(mfrow=c(2,1))
plot(routineTime[6000:8000],simulatedFlow[6000:8000,1],type="l")
points(eventData[[1]]$Time,eventData[[1]]$Flow,col="red")

plot(routineTime[6000:8000],simulatedTP[6000:8000,1],type="l")
points(eventData[[1]]$Time,eventData[[1]]$TP,col="red")


part=1
subpart=1
year=8
time<-matrix(seq(1/24,(365*20),1/24),ncol=20)
load(paste("~/Documents/code/Simulations/predicted/event/part",part_no,"/krigtp_part",part_no,"_subpart",subpart,".Rdata",sep=""))
event.krig.tp<-krig.tp
load(paste("~/Documents/code/Simulations/predicted/routine/part",part_no,"/krigtp_part",part_no,"_subpart",subpart,".Rdata",sep=""))
routine.krig.tp<-krig.tp
load(paste("~/Documents/code/Simulations/data/backtransformed_simulations/parts/tp/simulatedTP",part,".Rdata",sep=""))
real.<-matrix(simulatedTP[,subpart],ncol=20)
load(paste("data/simulated_sampled_data/event/eventDatapart",part,".Rdata",sep=""))



plot(routineTime[1:3000],simulatedTP[1:3000,72],type="l")
points(eventData[[72]]$Time,eventData[[72]]$TP,col="red")



plot(time[,year],event.krig.tp[[year]][[1]][],main=paste("Realisation number: ",subpart,". Year: ",year,sep=""),type="l",xlab="hours of year",ylab="TP",col="gray90",lty=4)
lines(time[,year],routine.krig.tp[[year]][[1]][],col="gray70",lty=4)
lines(time,simulatedTP[,subpart],col="black")
points(eventData[[subpart]]$Time,eventData[[subpart]]$TP,col="red")



plot(time[,2],real.[,2],type="l")
points(eventData[[1]]$Time,eventData[[1]]$TP,col="red")
lines(time[,2],event.krig.tp[[2]][[1]],col="red",lty=4)

plot(time[1:2000],real.[1:2000],type="l")
points(eventData[[1]]$Time,eventData[[1]]$TP,col="red")


part=round(no./250)+1
subpart=no.%%250
year<-sample(1:20,1)
load(paste("~/Documents/code/Simulations/predicted/event/part",part_no,"/krigtp_part",part_no,"_subpart",subpart,".Rdata",sep=""))
event.krig.tp<-krig.tp
load(paste("~/Documents/code/Simulations/predicted/routine/part",part_no,"/krigtp_part",part_no,"_subpart",subpart,".Rdata",sep=""))
routine.krig.tp<-krig.tp
load(paste("~/Documents/code/Simulations/data/backtransformed_simulations/parts/tp/simulatedTP",part,".Rdata",sep=""))
real.<-matrix(simulatedTP[,subpart],ncol=20)
plot(event.krig.tp[[year]][[1]][],main=paste("Realisation number: ",subpart,". Year: ",year,sep=""),type="l",xlab="hours of year",ylab="TP",col="blue",lty=4)
lines(routine.krig.tp[[year]][[1]][],col="red",lty=4)
lines(real.[,year],col="black")



load(paste("~/Documents/code/Simulations/predicted/routine/part1/krigtp_part1_subpart1.Rdata",sep=""))

###validation of variograms
library(gstat)
load("~/Documents/code/Simulations/temp_results/Verify_validation.Rdata")
g<-do.call(cbind,small.gstat)
fake<-small.gstat[[1]]
fake$gamma<-as.numeric(apply(g[,seq(3,length(names(g)),6)],1,mean))
plot(fake)
