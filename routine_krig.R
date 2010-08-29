checkout<-function(){
	ip=system("ifconfig -a | awk \'/eth0/{p=1}p&&/inet addr/{sub(\".*:\", \"\", $2);print $2;exit}\'",intern=T)

	load("routine_krigDB.Rdata")
	for(i in 1:(length(status[,1])+1)){
		if(i==(length(status[,1])+1))stop(NULL)
		if(is.na(status[i,3])||as.factor(ip)==status[i,3]){#The or part makes it check for previous work that wasn't finished.
			status[i,3]=ip
			save(status,file="routine_krigDB.Rdata")
			return(status[i,1:2])
		}
	}
}


checkin<-function(it.){
	load("routine_krigDB.Rdata")
	status[as.numeric(row.names(it.)),3]="done"
	save(status,file="routine_krigDB.Rdata")
}

###################################
###################################
##Before the real action starts lets make the routine time matrix
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
###Load the results of all the models
load("routine_based_model_summary.Rdata")

while(!is.null(part[,1])){
	part<-checkout()
	
	###load the sampled data - if needed
	if(temp!=part[,1]){
		load(paste("data/simulated_sampled_data/routine/routineFlowpart",part[,1],".Rdata",sep=""))#routineTP 250 realisations
		load(paste("data/simulated_sampled_data/routine/routineTPpart",part[,1],".Rdata",sep=""))#routineFlow 250 realisations
	}
	temp<-part[,1]
	###load the routinelikfit model
	load(paste("models/routine/part",part[,1],"/routinelikfit_part_",part[,1],"_subpart_",part[,2],".Rdata",sep=""))	
	###Load the simulated Flow for this realisation
	load(paste("data/backtransformed_simulations/parts/flow/ind/simFlow_part",part[,1],"_subpart",part[,2],".Rdata",sep=""))
	###Create a matrix of the flow, with each of the 20 columns represntting each year.
	real.flow<-matrix(simFlow,ncol=20)


###Get the lambda value of flow from the samples.
	flowlambda<-boxcox.fit(routineFlow[,part[,2]])$lambda
###Create a data frame of the sampled data of the appropriate realisation, with flow transformed (event has the lambda in the likfit model.
	obsdata<-data.frame(X=routineTime,Y=1,tp=routineTP[,part[,2]],flow=BCtransform(routineFlow[,part[,2]],flowlambda)$data)

###Now check if it will be a linear model or a lmcr. using the data.frame (results column spatial.mod)

	model<-(part[,1]-1)*250+part[,2] #convert the subpart and part numbers to a value between 1:2500

	if(results$spatial.mod[model]==1){ #Do the next step which uses the spatial model.
		###Create the geo object.
		obsdata.geo<-as.geodata(obsdata,coords.col=1:2,data.col=3,covar.col=4)
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
			###set the krigging model. - ordinary kriging, with the trend ~log(flow+0.001), using the likfit model of the realisation.
			krige.control<-krige.control(type.krige="SK",trend.d=trend.d,trend.l=trend.l,obj.model=routinelikfit[[1]])
			### Predict tp using the above model and save it to a list.
			krig.tp[[year]]<-krige.conv(obsdata.geo,locations=flow.[,1:2],krige=krige.control)
			
		}
		save(krig.tp,file=paste("predicted/routine/part",part[,1],"/krigtp_",part[,1],"_subpart",part[,2],".Rdata",sep=""))
		checkin(part)
	
	} else{ #Do the none spatial model (i.e do use a linear model to predict with.

		#Get the lambda value for tp from the routine likfit model
		tplambda<-summary(routinelikfit[[1]])$transformation[1,2]
		#transform the tp observations.	
		obsdata$tp<-BCtransform(obsdata$tp,tplambda)$data
		#fit a linear model to the data
		lin.model<-lm(tp~flow,data=obsdata)
		
		####Now do the predictions
		
		for(year in 1:20){
			#Make the predictions for the year...
			temp.predict<-predict.lm(lin.model,data.frame(flow=BCtransform(real.flow[,year],flowlambda)$data),se.fit=T)
			#Back transform the predictions
			krig.tp[[year]]<-backtransform.moments(tplambda,as.numeric(temp.predict$fit),as.numeric(temp.predict$se.fit))
		}
		save(krig.tp,file=paste("predicted/routine/part",part[,1],"/krigtp_",part[,1],"_subpart",part[,2],".Rdata",sep=""))
		checkin(part)	
	}
}

	


