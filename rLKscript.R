####################this script is for client based computers


####Get the distances
library(geoR)
date<-rep(seq.Date(as.Date("1991-01-01"),as.Date("1991-12-31"),by="days"),each=24)
#Extract the numeric day of each date
date<-as.numeric(format(date,"%d"))
#Join the day number with the hour of the day   
date<-paste(date,rep(seq(0:23),365)-1,sep=".")
#Then repeat the sequence above (hourly seq of one year) for 20 years
date<-rep(date,20)
###Create a custom time to represtent when the routine sample was taken (in days at hourly basis for 20 years)
routineTime<-(1:length(date))/24
routineTime<-routineTime[date==15.12]

###Set the check in and check out functions.
checkout<-function(){
	ip=system("ifconfig -a | awk \'/eth0/{p=1}p&&/inet addr/{sub(\".*:\", \"\", $2);print $2;exit}\'",intern=T)

	load("routineDB.Rdata")
	for(i in 1:(length(status[,1])+1)){
		if(i==(length(status[,1])+1))stop(NULL)
		if(is.na(status[i,3])||as.factor(ip)==status[i,3]){#The or part makes it check for previous work that wasn't finished.
			status[i,3]=ip
			save(status,file="routineDB.Rdata")
			return(status[i,1:2])
		}
	}
}


checkin<-function(it.){
	load("routineDB.Rdata")
	status[as.numeric(row.names(it.)),3]="done"
	save(status,file="routineDB.Rdata")
}


ini.=matrix(NA,ncol=2,nrow=6)
ini.[,1]=c(0.01,0.01,0.01,0.01,0.01,0.01)
ini.[,2]=c(100,20,30,10,150,300)
routinelikfit<-list()
#Have to intialise part for the while loop.
part<-data.frame(temp=1,temp=1)
##Create an indicator to see if the data set needs to be reloaded.
temp<-0

while(!is.null(part[,1])){
	###Checkout which one to do.	
	part<-checkout()
	
	###load the sampled data - if needed
	if(temp!=part[,1]){
		load(paste("data/simulated_sampled_data/routine/routineFlowpart",part[,1],".Rdata",sep=""))#routineTP 250 realisations
		load(paste("data/simulated_sampled_data/routine/routineTPpart",part[,1],".Rdata",sep=""))#routineFlow 250 realisations
	}
	temp<-part[,1]

	lambda<-boxcox.fit(routineTP[,part[,2]])$lambda
	flowlambda<-boxcox.fit(routineFlow[,part[,2]])$lambda

	routineOKData<-data.frame(x=routineTime,y=1,tp=routineTP[,part[,2]],flow=BCtransform(routineFlow[,part[,2]],flowlambda)$data)
	routineGeo<-as.geodata(routineOKData,coords.col=1:2,data.col=3,covar.col=4)

	routinelikfit[[1]]<-likfit(routineGeo,ini=ini.,fix.nugget = F,fix.lambda=T, lik.method = "REML",lambda=lambda,trend=trend.spatial(~flow))
	save(routinelikfit,file=paste("models/routine/part",part[,1],"/routinelikfit_part_",part[,1],"_subpart_",part[,2],".Rdata",sep=""))

	checkin(part)
}
