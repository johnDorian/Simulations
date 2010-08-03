setwd("~/Documents/code/Simulations/")
###Load the flow data again.
load("simflow5.Rdata")
load("eventlikfit_part5.Rdata")
load("eventData.Rdata")
eventData<-eventData[1001:1250]
krig.tp<-list()
library(geoR)

###Time
time<-matrix(seq(1/24,(365*20),1/24),ncol=20)
k=1
for(i in 1:250){ #i is for each realisation
### Seperate each realisation into individual years
	realisation<-matrix(temp[,i],ncol=20)

#######Get the obsereved samples into a geo object
###Get the lambda values
	flowlambda<-boxcox.fit(eventData[[i]]$Flow)$lambda
###Make the observation into a data frame with transformed values - dont transform TP as it is within the eventlikfit object.
	obdata<-data.frame(X=eventData[[i]]$Time,Y=1,tp=eventData[[i]]$TP,flow=BCtransform(eventData[[i]]$Flow,flowlambda)$data)
###Create the geo object.
	obdata<-obdata[!duplicated(obdata),]
	temp2<-as.geodata(obdata,coords.col=1:2,data.col=3,covar.col=4)
###State the trend of the krig model.
	trend.d<-trend.spatial(~flow,temp2)

	for(j in 1:20){	#j is for each year within each realisation.


###Get the data for the each year of the each realisation
		flow.<-data.frame(X=time[,j],Y=1,flw=BCtransform(realisation[,j],flowlambda)$data)
##Make this year a geodata object
		flow.sp<-as.geodata(flow.,coords.col=1:2,covar.col=3)
###Specify the krigging trend , based on all flow for the year.
		trend.l<-trend.spatial(~flw,flow.sp)


###set the krigging model. - ordinary kriging, with the trend ~log(flow+0.001), using the likfit model of the realisation.
		krige.control<-krige.control(type.krige="SK",trend.d=trend.d,trend.l=trend.l,obj.model=eventlikfit[[i]])
### Predict tp using the above model and save it to a list.
		krig.tp[[k]]<-krige.conv(temp2,locations=flow.[,1:2],krige=krige.control)
		k=k+1
	}
	if(i%%10==0){
		save(krig.tp,file=paste("krigtp_part5_subpart",i,".Rdata",sep=""))
		rm(krig.tp)
		gc()
		krig.tp<-list()
		k=1
	}
	print(i)
}


