####This first section is designed to automatically get the computer rolling.
###Assign a part and sub part number for this computer based on the database DB.Rdata - the script loops through the database gives its ip when checking a file out, then classifies the file as completed at the end.

##First section adds each computer to the client list.
#clients<-read.table("clients.txt",header=T)
#this.client<-system("ifconfig -a | awk \'/eth0/{p=1}p&&/inet addr/{sub(\".*:\", \"\", $2);print $2;exit}\'",intern=T)
#write.table(this.client,row.names=F,col.names=F,file="clients.txt",append=T)


###THe error checking of the checkout function is temporally changed to 1251 as we are only doing the second half of the dataset.
checkout<-function(){
	ip=system("ifconfig -a | awk \'/eth0/{p=1}p&&/inet addr/{sub(\".*:\", \"\", $2);print $2;exit}\'",intern=T)

	load("DB.Rdata")
	for(i in 1:(length(status[,1])+1)){
		if(i==(length(status[,1])+1))stop(NULL)
		if(is.na(status[i,3])||as.factor(ip)==status[i,3]){#The or part makes it check for previous work that wasn't finished.
			status[i,3]=ip
			save(status,file="DB.Rdata")
			return(status[i,1:2])
		}
	}
}


checkin<-function(it.){
	load("DB.Rdata")
	status[as.numeric(row.names(it.)),3]="done"
	save(status,file="DB.Rdata")
}


###Set up the intial values for the likfitting function
ini.<-matrix(NA,ncol=2,nrow=7)
ini.[,1]=c(0.1,0.1,0.1,0.1,0.1,0.1,0.1)
ini.[,2]=c(2,3,4,5,6,7,8)
###load the geoR library
library(geoR)
###Create a list for the results
eventlikfit<-list()
#Have to intialise part for the while loop.
part<-data.frame(temp=1,temp=1)
##Create an indicator to see if the data set needs to be reloaded.
temp<-0
while(!is.null(part[,1])){
	###Checkout which one to do.	
	part<-checkout()
	
	###load the sampled data - if needed
	if(temp!=part[,1]){
		load(paste("data/simulated_sampled_data/event/eventDatapart",part[,1],".Rdata",sep=""))
	}
	temp<-part[,1]

	lambda<-boxcox.fit(eventData[[part[,2]]]$TP)$lambda
	flowlambda<-boxcox.fit(eventData[[part[,2]]]$Flow)$lambda
#eventData$flow<-BCtransform(eventData$flow,flowlambda)$data
##Create a geodata object of the data.
#BCtransform(eventData[[i]]$Flow,flowlambda)$data
	eventOKData<-data.frame(x=eventData[[part[,2]]]$Time,y=1,tp=eventData[[part[,2]]]$TP,flow=BCtransform(eventData[[part[,2]]]$Flow,flowlambda)$data)
	eventOKData<-eventOKData[!duplicated(eventOKData[,1]),]
	eventGeo<-as.geodata(eventOKData,coords.col=1:2,data.col=3,covar.col=4)
###Plot a variogram of the data, using lambda to transform the data
#plot(variog(eventGeo,lambda=lambda,max.dist=30,trend=~flow))
###Use likfit to optimise the variogram model.
	eventlikfit[[1]]<-likfit(eventGeo,ini=ini.,fix.nugget = F,fix.lambda=T, lik.method = "REML",lambda=lambda,trend=trend.spatial(~flow))
#lines(eventlikfit)
	save(eventlikfit,file=paste("models/event/part",part[,1],"/eventlikfit_part_",part[,1],"_subpart_",part[,2],".Rdata",sep=""))
	checkin(part)
}



