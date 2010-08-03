####calculate annual load - 20 per realisation
####Calculate relaistaion rmse - 10 per file
####Do this 10 times
setwd("~/Documents/code/Simulations")
load("simtp10.Rdata")
tp=temp
load("simflow10.Rdata")
flow=temp
rm(temp)
gc()

load<-tp*flow
rm(tp)
gc()
annual.load<-matrix(NA,ncol=250,nrow=20)
for(i in 1:250){
	real.load<-matrix(load[,i],ncol=20)
	real.load<-colSums(real.load)
	annual.load[,i]<-real.load
}

library(geoR)

pred.load<-matrix(NA,ncol=250,nrow=20)
final=0
	for(file in 1:25){
##load the kriged data for the next ten realisations
		load(paste("krigtp_part10_subpart",file*10,".Rdata",sep=""))
		krig.counter=0
		for(realisation in 1:10){
			final=final+1
			flow.year<-matrix(flow[,final],ncol=20)
			for(i in 1:20){
				krig.counter=krig.counter+1
###predict must be changed for what ever the predicted object name actually is.
				pred.load[i,final]<-sum(flow.year[,i]*krig.tp[[krig.counter]]$predict)
			}
		}
	}
				
###Now we have two matrixs with annual load sums, 250 cols (for each realisation) and 20 rows (for each year of each realisation).
rmse <- function(obs, pred) sqrt(mean((obs-pred)^2))
temp<-NA
for(i in 1:250){
temp<-rbind(temp,rmse(annual.load[,i],pred.load[,i]))
}
temp<-temp[-1,]
hist(temp)

eventFinal<-list()
eventFinal[["ActualAnnualLoad"]]<-annual.load
eventFinal[["PredictedAnnualLoad"]]<-pred.load
eventFinal[["RMSE"]]<-as.numeric(temp)

save(eventFinal,file="eventfinal_part10.Rdata")
#load("eventfinal_part10.Rdata")






#################NOW to recreate above but for every part and subpart of the realisations for the event based load estimation. - experimental as of 3/8/10. 



#setwd("~/Documents/code/Simulations")

for(part_no in 1:10){

	load(paste("~/Documents/code/Simulations/simtp",part_no,".Rdata",sep=""))
	tp=temp
	load(paste("~/Documents/code/Simulations/simflow",part_no,".Rdata",sep=""))
	flow=temp
	rm(temp)
	gc()

	load<-tp*flow
	rm(tp)
	gc()
	annual.load<-matrix(NA,ncol=250,nrow=20)
	for(i in 1:250){
		real.load<-matrix(load[,i],ncol=20)
		real.load<-colSums(real.load)
		annual.load[,i]<-real.load
	}

	library(geoR)

	pred.load<-matrix(NA,ncol=250,nrow=20)
	final=0
		for(file in 1:25){
##load the kriged data for the next ten realisations
			load(paste("~/Documents/code/Simulations/predicted/event/part",part_no,"/krigtp_part",part_no,"_subpart",file*10,".Rdata",sep=""))
			krig.counter=0
			for(realisation in 1:10){
				final=final+1
				flow.year<-matrix(flow[,final],ncol=20)
				for(i in 1:20){
					krig.counter=krig.counter+1
###predict must be changed for what ever the predicted object name actually is.
					pred.load[i,final]<-sum(flow.year[,i]*krig.tp[[krig.counter]]$predict)
				}
			}
		}
				
###Now we have two matrixs with annual load sums, 250 cols (for each realisation) and 20 rows (for each year of each realisation).
	rmse <- function(obs, pred) sqrt(mean((obs-pred)^2))
	temp<-NA
	for(i in 1:250){
		temp<-rbind(temp,rmse(annual.load[,i],pred.load[,i]))
	}
	temp<-temp[-1,]
	hist(temp)

	eventFinal<-list()
	eventFinal[["ActualAnnualLoad"]]<-annual.load
	eventFinal[["PredictedAnnualLoad"]]<-pred.load
	eventFinal[["RMSE"]]<-as.numeric(temp)

	save(eventFinal,file=paste("~/Documents/code/Simulations/eventfinal_part",part_no,".Rdata",sep=""))
}
