###An inestigation into the problem with realisation 1871 year 7 of the event sampled data.
realisation=1871
part=floor(realisation/250)
subpart=realisation%%250

###Load the actual data
load(paste('~/Documents/code/Simulations/data/backtransformed_simulations/parts/tp/simulatedTP',part,'.Rdata',sep=""))
load(paste('~/Documents/code/Simulations/data/backtransformed_simulations/parts/flow/simulatedFlow',part,'.Rdata',sep=""))
###Pull out the relevant realisation
simulatedTP<-simulatedTP[,subpart]
simulatedFlow<-simulatedFlow[,subpart]
###Load the predicted data
load(paste('~/Documents/code/Simulations/predicted/event/part',part,'/krigtp_part',part,'_subpart',subpart,'.Rdata',sep=""))
###Load the likfit model
load(paste('~/Documents/code/Simulations/models/event/part',part,"/eventlikfit_part_",part,'_subpart_',subpart,'.Rdata',sep=""))



###Re-run the predictions
###load the geoR library
library(geoR)
###Create a list for the results
krig.tp2<-list()
###Create a matrix of the time for each ob for the 20 years
time<-matrix(seq(1/24,(365*20),1/24),ncol=20)


load(paste("data/simulated_sampled_data/event/eventDatapart",part,".Rdata",sep=""))# 250 realisations of sampled data. - eventData[[1-250]]
###load the routinelikfit model
load(paste("models/event/part",part,"/eventlikfit_part_",part,"_subpart_",subpart,".Rdata",sep=""))#eventlikfit - eventlikfit[[1]] only.
###Load the simulated Flow for this realisation
load(paste("data/backtransformed_simulations/parts/flow/ind/simFlow_part",part,"_subpart",subpart,".Rdata",sep=""))
###Create a matrix of the flow, with each of the 20 columns represntting each year.
real.flow<-matrix(simFlow,ncol=20)


###Get the lambda value of flow from the samples.
flowlambda<-boxcox.fit(eventData[[subpart]]$Flow)$lambda
###Create a data frame of the sampled data of the appropriate realisation, with flow transformed (event has the lambda in the likfit model.
obsdata<-data.frame(X=eventData[[subpart]]$Time,Y=1,tp=eventData[[subpart]]$TP,flow=BCtransform(eventData[[subpart]]$Flow,flowlambda)$data)
###Get rid of duplicated observations. (These are due to the merging of routine and event samples)
obsdata<-obsdata[!duplicated(obsdata),]
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
	rm(flow.sp)
	gc()
	###set the krigging model. - ordinary kriging, with the trend ~log(flow+0.001), using the likfit model of the realisation.
	krige.control<-krige.control(type.krige="SK",trend.d=trend.d,trend.l=trend.l,obj.model=eventlikfit[[1]])
	### Predict tp using the above model and save it to a list.
	krig.tp2[[year]]<-krige.conv(obsdata.geo,locations=flow.[,1:2],krige=krige.control)
		
}
#save(krig.tp,file=paste("predicted/event/part",part[,1],"/krigtp_",part[,1],"_subpart",part[,2],".Rdata",sep=""))

#Check to see if there is a difference...
load(paste('~/Documents/code/Simulations/models/event/part',part,'/eventlikfit_part_',part,'_subpart_',subpart,'.Rdata',sep=''))
plot(krig.tp[[6]][[1]],krig.tp2[[6]][[1]])
