partno<-commandArgs(TRUE)
if(!is.numeric(partno))stop("Input must be numeric")
if(partno<0||partno>10)stop("Input must be between 1 and 10")



load(paste("~/Documents/code/Simulations/data/simulated_sampled_data/event/eventDatapart",partno,".Rdata"))
library(geoR)
eventlikfit<-list()

ini.<-matrix(NA,ncol=2,nrow=13)
ini.[,1]=c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)
ini.[,2]=c(2,3,4,5,6,7,8,10,20,30,50,70,100)

for(i in 1:250){


	lambda<-boxcox.fit(eventData[[i]]$TP)$lambda
	flowlambda<-boxcox.fit(eventData[[i]]$Flow)$lambda
#eventData$flow<-BCtransform(eventData$flow,flowlambda)$data
##Create a geodata object of the data.
#BCtransform(eventData[[i]]$Flow,flowlambda)$data
	eventOKData<-data.frame(x=eventData[[i]]$Time,y=1,tp=eventData[[i]]$TP,flow=BCtransform(eventData[[i]]$Flow,flowlambda)$data)
	eventOKData<-eventOKData[!duplicated(eventOKData),]
	eventGeo<-as.geodata(eventOKData,coords.col=1:2,data.col=3,covar.col=4)
###Plot a variogram of the data, using lambda to transform the data
#plot(variog(eventGeo,lambda=lambda,max.dist=30,trend=~flow))
###Use likfit to optimise the variogram model.
		eventlikfit[[i]]<-likfit(eventGeo,ini=ini.,fix.nugget = F,fix.lambda=T, lik.method = "REML",lambda=lambda,trend=trend.spatial(~flow))
#lines(eventlikfit)
print(i)
}
save(eventlikfit,file="eventlikfit_part10.Rdata")
rm(eventlikfit)
gc()

