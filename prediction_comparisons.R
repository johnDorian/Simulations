###Now  it's time to find out which one is the best
hours.above<-list()
rmse.results<-matrix(NA,ncol=2,nrow=2500)
rmse <- function(obs, pred) sqrt(mean((obs-pred)^2))
for(part in 1:10){
	load(paste("data/backtransformed_simulations/parts/tp/simulatedTP",part,".Rdata",sep=""))
	for(subpart in 1:250){
		real.tp<-matrix(simulatedTP[,subpart],ncol=20)
		load(paste("predicted/event/part",part,"/krigtp_part",part,"_subpart",subpart,".Rdata",sep=""))
		event.tp<-matrix(NA,ncol=20,nrow=8760)
		routine.tp<-matrix(NA,ncol=20,nrow=8760)
		for(i in 1:20){
			event.tp[,i]<-krig.tp[[i]]$predict
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



plot(real.tp[,2],type="l")
lines(event.tp[,2],col="red")
lines(routine.tp[,2],col="blue")


days=data.frame(actual=rep(NA,20),event=NA,routine=NA)
for(i in 1:20){
	days[i,1]<-length(real.tp[real.tp[,i]>0.2,1])
	days[i,2]<-length(event.tp[event.tp[,i]>0.2,1])
	days[i,3]<-length(routine.tp[routine.tp[,i]>0.2,1])
}


rmse.results[1,1:2]<-c(rmse(days[,1],days[,2]),rmse(days[,1],days[,3]))

pdf("RMSE_hist.pdf",encoding="UTF-8")
par(mfrow=c(2,1))
hist(rmse.results[,1]/24,main="RMSE of threshold exceedance",xlab="Event-based sampling")
text(5.2,600,paste("\u03bc =",round(mean(rmse.results[,1]/24)),",\u03c3 =",round(sd(rmse.results[,1]/24))))
hist(rmse.results[,2]/24,main="",xlab="Routine-based sampling")
text(10,400,paste("\u03bc =",round(mean(rmse.results[,2]/24)),",\u03c3 =",round(sd(rmse.results[,2]/24))))
dev.off()
pdf("RMSE_dots.pdf")
re<-rgb(red=255,blue=0,green=0,max=255,alpha=120)
bla<-rgb(red=0,blue=0,green=0,max=255,alpha=120)
ylim=c(min(rmse.results/24),max(rmse.results/24))
plot(rmse.results[,1]/24,ylab="RMSE exceedance",xlab="Realisation",ylim=ylim,pch=16,col=bla,main="RMSE of exceednace for each realisation of each samping sheme")
points(rmse.results[,2]/24,col=re,pch=16)
legend("topleft",c("Event","routine"),col=c(bla,re),pch=16,bg="white")
dev.off()



####Now get the mean amount of days above the threshold
##Convert the list into a dataframe
cb<-do.call(cbind,hours.above)
##Every 3rd column from the first column is the real values
re.seq<-seq(1,length(cb),3)
##Every 3rd column from the first column is the event values
eve.seq<-seq(2,length(cb),3)
##Every 3rd column from the first column is the routine values
rou.seq<-seq(3,length(cb),3)

##Get the average days above the threshold of the real data.
averageAbove<-data.frame(actual=as.numeric(colSums(cb[,re.seq])/20/24),event=as.numeric(colSums(cb[,eve.seq])/20/24),routine=as.numeric(colSums(cb[,rou.seq])/20/24))
##Now create a boxplot of the results
pdf("bplot_average_exc.pdf")
boxplot(averageAbove,main="Average days above threshold per year per realisation",xlab="Sample type",ylab="Days above 0.2")
dev.off()
###Save the data, so I don't have to do this again.
save(averageAbove,file="averageDaysAbove0.2.Rdata")
save(rmse.results,file="averageAbovehoursRMSE.Rdata")

########################Now lets do the load...



library (geoR)


for (part_no in 1:10){
	###Load the tp and flow data (complete simulated)
	load(paste("~/Documents/code/Simulations/data/backtransformed_simulations/parts/tp/simulatedTP",part_no,".Rdata",sep=""))

	load(paste("~/Documents/code/Simulations/data/backtransformed_simulations/parts/flow/simflow",part_no,".Rdata",sep=""))
	###Calulate the load
	load<-simulatedTP*simulatedFlow
	###remove the tp - don't need it anymore
	rm(simulatedTP)
	gc()
	###Cycle through the 250 realisations of the simulated data and get the annual load
	annual.load<-matrix(NA,ncol=250,nrow=20)
	for(i in 1:250){
		real.load<-matrix(load[,i],ncol=20)
		real.load<-colSums(real.load)
		annual.load[,i]<-real.load
	}

	
	###Create a matrix for predicted load 
	pred.load<-matrix(NA,ncol=250,nrow=20)
	final=0
		for (subpart in 1:250){
##load the kriged data for the next ten realisations
			load(paste("~/Documents/code/Simulations/predicted/event/part",part_no,"/krigtp_part",part_no,"_subpart",subpart,".Rdata",sep=""))
			krig.counter=0
			for (realisation in 1:10){
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











