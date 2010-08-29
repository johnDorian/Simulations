setwd("~/Documents/code/Simulations/models/event/")
library(geoR)
spatial<-vector()
ns.spatial<-vector()
range<-vector()
real<-vector()
nugget<-vector()
sill<-vector()
i=0
for(part in 1:10){
	for(subpart in 1:250){
		i=i+1
		load(paste("part",part,"/eventlikfit_part_",part,"_subpart_",subpart,".Rdata",sep=""))
		spatial[i]<-summary(eventlikfit[[1]])$likelihood$AIC
		ns.spatial[i]<-summary(eventlikfit[[1]])$ns.likelihood$AIC
		range[i]<-summary(eventlikfit[[1]])$spatial.component[2,2]
		sill[i]<-summary(eventlikfit[[1]])$spatial.component[1,2]
		nugget[i]<-summary(eventlikfit[[1]])$nugget.component[1,2]

		if(spatial[i]<ns.spatial[i])real[i]<-1 else real[i]<-0
	}
	print(part)
}

hist(range) # the range of each model.
sum(real) # If the spatial model is better then it is 1, therefore this indicates how many models are spatial.
results<-data.frame(realisation=seq(1,2500),spatial.AIC=spatial,ns.spatial.AIC=ns.spatial,range,sill,nugget,spatial.mod=real)
save(results,file="~/Documents/code/Simulations/event_based_model_summary.Rdata")

setwd("~/Documents/code/Simulations/models/routine/")
library(geoR)
spatial<-vector()
ns.spatial<-vector()
range<-vector()
real<-vector()
nugget<-vector()
sill<-vector()
i=0
for(part in 1:10){
	for(subpart in 1:250){
		i=i+1
		load(paste("part",part,"/routinelikfit_part_",part,"_subpart_",subpart,".Rdata",sep=""))
		spatial[i]<-summary(routinelikfit[[1]])$likelihood$AIC
		ns.spatial[i]<-summary(routinelikfit[[1]])$ns.likelihood$AIC
		range[i]<-summary(routinelikfit[[1]])$spatial.component[2,2]
		sill[i]<-summary(routinelikfit[[1]])$spatial.component[1,2]
		nugget[i]<-summary(routinelikfit[[1]])$nugget.component[1,2]
		#for this to be a spatial model it must also have a range less than 365 days.
		if(spatial[i]<ns.spatial[i]&&range[i]<365)real[i]<-1 else real[i]<-0
	}
	print(part)
}

hist(results[results[,"spatial.mod"]==1,"range"]) # the range of each spatial model.
sum(real) # If the spatial model is better then it is 1, therefore this indicates how many models are spatial.
results<-data.frame(realisation=seq(1,2500),spatial.AIC=spatial,ns.spatial.AIC=ns.spatial,range,sill,nugget,spatial.mod=real)
save(results,file="~/Documents/code/Simulations/routine_based_model_summary.Rdata")



