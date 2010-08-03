setwd("~/Documents/code/Simulations/")
###Load the flow data again.
load("simflow1.Rdata")
#TODO:Make simflow2 work some how. #load("simflow2.Rdata") if(i==500)load(sim2)
load("routinelikfit1.Rdata")
load("routineLFsummary.Rdata")
load("routineFlow.Rdata")
load("routineTP.Rdata")



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


krig.tp<-list()
library(geoR)

###Time
time<-matrix(seq(1/24,(365*20),1/24),ncol=20)
k=1
for(i in 1:2){ #i is for each realisation
temp.it<-NA
if(i>250&&i<500)temp.it=i%%250 else if(i%%250==0) temp.it=250 else if(i<250)temp.it=i


###---------- Load the second file for the simflow data---------#####
	if(i==251){
		rm(temp)
		gc()		
		load("simflow2.Rdata")
	}
###--------------------------------------------------------------####

	if(routine_likfit_summary$Spatial_binomial[[i]]==1){
### Seperate each realisation into individual years

		realisation<-matrix(temp[,temp.it],ncol=20)

#######Get the obsereved samples into a geo object
###Get the lambda values
		flowlambda<-boxcox.fit(routineFlow[,i])$lambda
###Make the observation into a data frame with transformed values - dont transform TP as it is within the eventlikfit object.
		obdata<-data.frame(X=routineTime,Y=1,tp=routineTP[,i],flow=BCtransform(routineFlow[,i],flowlambda)$data)
###Create the geo object.
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
			krige.control<-krige.control(type.krige="SK",trend.d=trend.d,trend.l=trend.l,obj.model=routinelikfit[[i]])
### Predict tp using the above model and save it to a list.
			krig.tp[[k]]<-krige.conv(temp2,locations=flow.[,1:2],krige=krige.control)
			k=k+1
		}
	}else{
###Get the lambda values
		flowlambda<-boxcox.fit(routineFlow[,i])$lambda
		tplambda<-boxcox.fit(routineTP[,i])$lambda
###Make the observation into a data frame with transformed values - dont transform TP as it is within the eventlikfit object.
		obdata<-data.frame(X=routineTime,Y=1,tp=BCtransform(routineTP[,i],tplambda)$data,flow=BCtransform(routineFlow[,i],flowlambda)$data)
		lin.model<-lm(tp~flow,data=obdata)
###Make a prediction 'grid' of flow

		temp.flow<-matrix(BCtransform(temp[,temp.it],flowlambda)$data,ncol=20)
		
###loop through each year - not nessacry, but makes things easier if each method predicts on a an annual basis
		for(j in 1:20){
###Create a dataframe for each year within each realisation
			pred.mat<-data.frame(flow=temp.flow[,j])		
###Make some predictions
			temp.predict<-predict(lin.model,pred.mat,se.fit=T)
###Back transform the preditions and save them to the list
			krig.tp[[k]]<-backtransform.moments(lambda=tplambda,as.numeric(temp.predict$fit),as.numeric(temp.predict$se.fit))
			k=k+1
		}
	}



	if(i%%10==0){
		save(krig.tp,file=paste("routine_predicted_1_subpart",i,".Rdata",sep=""))
		rm(krig.tp)
		gc()
		krig.tp<-list()
		k=1
	}
	print(i)
}



##
