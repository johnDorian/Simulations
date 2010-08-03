####calculate annual load - 20 per realisation
####Calculate relaistaion rmse - 10 per file
####Do this 10 times
library(geoR)
for(part_no in 1:10){

	load(paste("~/Documents/code/Simulations/data/backtransformed_simulations/parts/tp/simtp",part_no,".Rdata",sep=""))
	tp=temp
	load(paste("~/Documents/code/Simulations/data/backtransformed_simulations/parts/flow/simflow",part_no,".Rdata",sep=""))
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
#	hist(temp)

	eventFinal<-list()
	eventFinal[["ActualAnnualLoad"]]<-annual.load
	eventFinal[["PredictedAnnualLoad"]]<-pred.load
	eventFinal[["RMSE"]]<-as.numeric(temp)

	save(eventFinal,file=paste("~/Documents/code/Simulations/temp_results/eventfinal_part",part_no,".Rdata",sep=""))
}


####And now for the routine dataset.

for(part_no in 1:5){
adj_part_no=0
	for(rep in 1:2){##This for loop is so the 250 realisations of tp and flow are combined to make 500 realisation objects like the same as the length of the routine objects.
		adj_part_no=adj_part_no+1
		if(adj_part_no==1){
			load(paste("~/Documents/code/Simulations/data/backtransformed_simulations/parts/tp/simtp",adj_part_no,".Rdata",sep=""))
			tp=temp
			load(paste("~/Documents/code/Simulations/data/backtransformed_simulations/parts/flow/simflow",adj_part_no,".Rdata",sep=""))
			flow=temp
			rm(temp)
			gc()
		}else{
			load(paste("~/Documents/code/Simulations/data/backtransformed_simulations/parts/tp/simtp",adj_part_no,".Rdata",sep=""))
			tp=cbind(tp,temp)
			load(paste("~/Documents/code/Simulations/data/backtransformed_simulations/parts/flow/simflow",adj_part_no,".Rdata",sep=""))
			flow=cbind(flow,temp)
			rm(temp)
			gc()
		}
	}

	load<-tp*flow
	rm(tp)
	gc()
	annual.load<-matrix(NA,ncol=500,nrow=20)
	for(i in 1:500){
		real.load<-matrix(load[,i],ncol=20)
		real.load<-colSums(real.load)
		annual.load[,i]<-real.load
	}

	

	pred.load<-matrix(NA,ncol=500,nrow=20)
	final=0
		for(file in 1:50){
##load the kriged data for the next ten realisations
			load(paste("~/Documents/code/Simulations/predicted/routine/part",part_no,"/routine_predicted_",part_no,"_subpart",file*10,".Rdata",sep=""))
			krig.counter=0
			for(realisation in 1:10){
				final=final+1
				flow.year<-matrix(flow[,final],ncol=20)
				for(i in 1:20){
					krig.counter=krig.counter+1
###Need to change this section as there is both lists and kriged
					if(class(krig.tp[[krig.counter]])=="kriging"){
						pred.load[i,final]<-sum(flow.year[,i]*krig.tp[[krig.counter]]$predict)
					}else{
						pred.load[i,final]<-sum(flow.year[,i]*krig.tp[[krig.counter]]$mean)
					}
				}
			}
		}
				
###Now we have two matrixs with annual load sums, 250 cols (for each realisation) and 20 rows (for each year of each realisation).
	rmse <- function(obs, pred) sqrt(mean((obs-pred)^2))
	temp<-NA
	for(i in 1:500){
		temp<-rbind(temp,rmse(annual.load[,i],pred.load[,i]))
	}
	temp<-temp[-1,]
#	hist(temp)

	routineFinal<-list()
	routineFinal[["ActualAnnualLoad"]]<-annual.load
	routineFinal[["PredictedAnnualLoad"]]<-pred.load
	routineFinal[["RMSE"]]<-as.numeric(temp)

	save(routineFinal,file=paste("~/Documents/code/Simulations/temp_results/routinefinal_part",part_no,".Rdata",sep=""))
}




####Now i have to add an extra bit of code to join them al up into two files. one for routine and one for event.
##cbind works for the first two, but not for the last one - use matrix to fix it up.
load(paste("~/Documents/code/Simulations/temp_results/eventfinal_part1.Rdata",sep=""))
temp<-eventFinal
for(part_no in 2:10){
	load(paste("~/Documents/code/Simulations/temp_results/eventfinal_part",part_no,".Rdata",sep=""))
	temp[[1]]<-cbind(temp[[1]],eventFinal[[1]])
	temp[[2]]<-cbind(temp[[2]],eventFinal[[2]])
	temp[[3]]<-cbind(temp[[3]],eventFinal[[3]])
}
temp[[3]]<-matrix(temp[[3]],ncol=1)
eventFinalAll<-temp
save(eventFinalAll,file="~/Documents/code/Simulations/temp_results/eventFinalAll.Rdata")


##Routine
load(paste("~/Documents/code/Simulations/temp_results/routinefinal_part1.Rdata",sep=""))
temp<-routineFinal
for(part_no in 2:5){
	load(paste("~/Documents/code/Simulations/temp_results/routinefinal_part",part_no,".Rdata",sep=""))
	temp[[1]]<-cbind(temp[[1]],routineFinal[[1]])
	temp[[2]]<-cbind(temp[[2]],routineFinal[[2]])
	temp[[3]]<-cbind(temp[[3]],routineFinal[[3]])
}
temp[[3]]<-matrix(temp[[3]],ncol=1)
routineFinalAll<-temp
save(routineFinalAll,file="~/Documents/code/Simulations/temp_results/routineFinalAll.Rdata")



###Just for kicks
load("~/Documents/code/Simulations/temp_results/eventFinalAll.Rdata")
load("~/Documents/code/Simulations/temp_results/routineFinalAll.Rdata")
hist(eventFinalAll[[3]]/1000)
dev.new()
hist(routineFinalAll[[3]]/1000)

summary(eventFinalAll[[3]]/1000)
summary(routineFinalAll[[3]]/1000)

