checkout<-function(){
	ip=system("ifconfig -a | awk \'/eth0/{p=1}p&&/inet addr/{sub(\".*:\", \"\", $2);print $2;exit}\'",intern=T)

	load("event_krigDB.Rdata")
	for(i in 1:(length(status[,1])+1)){
		if(i==(length(status[,1])+1))stop(NULL)
		if(is.na(status[i,3])||as.factor(ip)==status[i,3]){#The or part makes it check for previous work that wasn't finished.
			status[i,3]=ip
			save(status,file="event_krigDB.Rdata")
			return(status[i,1:2])
		}
	}
}


checkin<-function(it.){
	load("event_krigDB.Rdata")
	status[as.numeric(row.names(it.)),3]="done"
	save(status,file="event_krigDB.Rdata")
}

###################################
###################################
###load the geoR library
library(geoR)
###Create a list for the results
krig.tp<-list()
#Have to intialise part for the while loop.
part<-data.frame(temp=1,temp=1)
##Create an indicator to see if the data set needs to be reloaded.
temp<-0
###Create a matrix of the time for each ob for the 20 years
time<-matrix(seq(1/24,(365*20),1/24),ncol=20)


while(!is.null(part[,1])){
	part<-checkout()
	
	###load the sampled data - if needed
	if(temp!=part[,1]){
		load(paste("data/simulated_sampled_data/event/eventDatapart",part[,1],".Rdata",sep=""))# 250 realisations of sampled data. - eventData[[1-250]]
	}
	temp<-part[,1]
	###load the routinelikfit model
	load(paste("models/event/part",part[,1],"/eventlikfit_part_",part[,1],"_subpart_",part[,2],".Rdata",sep=""))#eventlikfit - eventlikfit[[1]] only.
	###Load the simulated Flow for this realisation
	load(paste("data/backtransformed_simulations/parts/flow/ind/simFlow_part",part[,1],"_subpart",part[,2],".Rdata",sep=""))
	###Create a matrix of the flow, with each of the 20 columns represntting each year.
	real.flow<-matrix(simFlow,ncol=20)
	rm(simFlow)
	gc()


###Get the lambda value of flow from the samples.
	flowlambda<-boxcox.fit(eventData[[part[,2]]]$Flow)$lambda
###Create a data frame of the sampled data of the appropriate realisation, with flow transformed (event has the lambda in the likfit model.
	obsdata<-data.frame(X=eventData[[part[,2]]]$Time,Y=1,tp=eventData[[part[,2]]]$TP,flow=BCtransform(eventData[[part[,2]]]$Flow,flowlambda)$data)
###Get rid of duplicated observations. (These are due to the merging of routine and event samples)
	obsdata<-obsdata[!duplicated(obsdata),]
	gc()

###Create the geo object.
	obsdata.geo<-as.geodata(obsdata,coords.col=1:2,data.col=3,covar.col=4)
	rm(obsdata)
	gc()
	###State the trend of the krig model.
	trend.d<-trend.spatial(~flow,obsdata.geo)

	###############Now do the predictions using 
	for(year in 1:20){	#j is for each year within each realisation.
	
		###Get the data for the each year of the each realisation
		flow.<-data.frame(X=time[,year],Y=1,flw=BCtransform(real.flow[,year],flowlambda)$data)
		##Make this year a geodata object
		flow.sp<-as.geodata(flow.,coords.col=1:2,covar.col=3)

		###Specify the krigging trend , based on all flow for the year.
		trend.l<-trend.spatial(~flw,flow.sp)
		rm(flow.sp)
		gc()
		###set the krigging model. - ordinary kriging, with the trend ~log(flow+0.001), using the likfit model of the realisation.
		krige.control<-krige.control(type.krige="SK",trend.d=trend.d,trend.l=trend.l,obj.model=eventlikfit[[1]])
		### Predict tp using the above model and save it to a list.
		krig.tp[[year]]<-krige.conv(obsdata.geo,locations=flow.[,1:2],krige=krige.control)
		
	}
	save(krig.tp,file=paste("predicted/event/part",part[,1],"/krigtp_",part[,1],"_subpart",part[,2],".Rdata",sep=""))
	checkin(part)

}
