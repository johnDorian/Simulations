setwd("~/Documents/code/Simulations")
###Create a seq to represent the realisations
seq<-1:2500
#how many samples to take from the population
n<-250
###Sample the realisation
random.cols<-sample(x=seq,size=n,replace=FALSE)
##Create a seq. to seperate the samples for each file
breaks<-seq(0,2500,by=250)
##Create a list and seperate each file realisations out. and modify the value so it represents the subpart of the file.
real.<-list()
for(i in 1:(length(breaks)-1)){
	real.[[i]]<-sort(random.cols[random.cols<breaks[i+1]& random.cols>breaks[i]]%%250)
}

### Load the gstat library
library(gstat);library(geoR)
sim=24*365*20
date=(seq(1,sim)-1)*0.04166667
small.gstat<-list()
	
output=1

for(i in 1:10){
	load(paste("data/simulated_raw/simtp/simtp",i,".Rdata",sep=""))#sim.tp
	load(paste("data/simulated_raw/simflow/simflow",i,".Rdata",sep=""))#sim.flow

	for(it in 1:length(real.[[i]])){
		data <- data.frame(X=date,Y=1,TP=sim.tp[,real.[[i]]][it],FLOW=sim.flow[,real.[[i]]][it])
### Create gstat object of the data
		spdf <- SpatialPointsDataFrame(data[,1:2],data)
### Create A gstat object with TP and flow.
		g = gstat(NULL, "TP", TP ~ 1,spdf,maxdist=200,dummy=TRUE,beta=mean(data$TP))
		g = gstat(g, "FLOW", FLOW ~ 1,spdf,maxdist=200,dummy=TRUE,beta=mean(data$FLOW))
		small.gstat[[i]] = variogram(g,cutoff=150,width=10)
		print(output)
		output=output+1
	}
}

save(small.gstat,file="sim_subset_gstat_validation.Rdata")
###Now to have a look at the results
