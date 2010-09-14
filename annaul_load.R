


realisation=0 #An iterator for 1:2500 realisations
###Now some matrixs for the results - 2500 realisations, and 20 years in each realisation
eventTP<-matrix(nrow=20,ncol=2500)
routineTP<-matrix(nrow=20,ncol=2500)
actualTP<-matrix(nrow=20,ncol=2500)
actualFlow<-matrix(nrow=20,ncol=2500)
#Firstly we will get the annual load for the predicted things
for(part in 1:10){
#realisation=realisation+1
load(paste('~/Documents/code/Simulations/data/backtransformed_simulations/parts/tp/simulatedTP',part,'.Rdata',sep=""))
load(paste('~/Documents/code/Simulations/data/backtransformed_simulations/parts/flow/simulatedFlow',part,'.Rdata',sep=""))
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
		##Go through this relisation and get the sum of each year.
		for(year in 1:20){
			eventTP[year,realisation]<-sum(eventkrig[[year]][[1]])
			routineTP[year,realisation]<-sum(krig.tp[[year]][[1]])
		}

	}
print(realisation)
}
##Now calculate the actual loads...
actualLoad<-actualTP*actualFlow
eventLoad<-eventTP*actualFlow
routineLoad<-routineTP*actualFlow





###############################
#####This next section will look at how many hours in the year the tp is over the threshold of 0.2
###############################
for(part in 1:10){
	load(paste("data/backtransformed_simulations/parts/tp/simulatedTP",part,".Rdata",sep=""))
	for(subpart in 1:250){
		real.tp<-matrix(simulatedTP[,subpart],ncol=20)


###Old stuff below
hours.above<-list()
rmse.results<-matrix(NA,ncol=2,nrow=2500)
rmse <- function(obs, pred) sqrt(mean((obs-pred)^2))
for(part in 1:10){
	load(paste("data/backtransformed_simulations/parts/tp/simulatedTP",part,".Rdata",sep=""))
	for(subpart in 1:250){
		real.tp<-matrix(simulatedTP[,subpart],ncol=20)
		load(paste("predicted/event/part",part,"/krigtp_part",part,"_subpart",subpart,".Rdata",sep=""))
		event.tp<-matrix(NA,ncol=20,nrow=8760)
		for(i in 1:20){
			event.tp[,i]<-krig.tp[[i]][[1]]
		}
		routine.tp<-matrix(NA,ncol=20,nrow=8760)
		load(paste("predicted/routine/part",part,"/krigtp_",part,"_subpart",subpart,".Rdata",sep=""))		
		for(i in 1:20){
			routine.tp[,i]<-krig.tp[[i]][[1]]
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




