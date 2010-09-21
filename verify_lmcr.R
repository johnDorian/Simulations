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
	real.[[i]]<-sort(random.cols[random.cols<=breaks[i+1]& random.cols>breaks[i]]%%250)
}

###Make sure the the amount of realisations to do is the same as what was requested.
if (sum(do.call(cbind,lapply(real.,length)))!=n) stop("Error in code")

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
		data <- data.frame(X=date,Y=1,TP=sim.tp[,real.[[i]][it]],FLOW=sim.flow[,real.[[i]][it]])
### Create gstat object of the data
		spdf <- SpatialPointsDataFrame(data[,1:2],data)
### Create A gstat object with TP and flow.
		g = gstat(NULL, "TP", TP ~ 1,spdf,maxdist=200,dummy=TRUE,beta=mean(data$TP))
		g = gstat(g, "FLOW", FLOW ~ 1,spdf,maxdist=200,dummy=TRUE,beta=mean(data$FLOW))
		small.gstat[[output]] = variogram(g,cutoff=150,width=10)
		print(output)
		output=output+1
	}
	small.gstat[[output]]<-real.
	save(small.gstat,file="sim_subset_gstat_validation.Rdata")
}

###Now to have a look at the results

###Now let's look at the results of the variograms.
#TODO: fix the next bit up. Make a nice plotting fuinction for lmcr, for verification and the actual modeled data.

##load the lmcrs of the simulated data
load('~/Documents/code/Simulations/sim_subset_gstat_validation.Rdata')
#Load the 
load('~/Documents/code/Simulations/data/lmcr_model/lmcr.Rdata')
############################################
##Reformat the data for plotting purposes.
############################################
##Get the distance values for plotting
dist=small.gstat[[1]]$dist[1:(length(small.gstat[[1]]$dist)/3)]
##Change the list into a matrix.
test1<-do.call(cbind,small.gstat[-248])
##Get the gamma values from the list
gamma<-test1[,seq(3,length(test1),by=6)]
##Get the mean of the gamma values
mean.gamma<-apply(gamma,1,mean)
##Get the standard deviation from the gamma values.
sd.gamma<-apply(gamma,1,sd)
###Get the length of simulations that where sampled.
n<-length(gamma[1,])
##Find the error of each bin.
error <- qnorm(0.975)*sd.gamma/sqrt(n)
##Get the lower and upper CI for each bin.
lower <- mean.gamma-error
upper <- mean.gamma+error


########################################
###Plotting time. ######################
########################################
###Below is the make up for my on function to plot lmcr's. 


##Now put the CI and the mean values into a matrix for each variogram type order=(TP.FLOW, FLOW,TP) - this is if the gstat tricking dosent work.
mean.gamma<-matrix(mean.gamma,nrow=15)
lower<-matrix(lower,nrow=15)
upper<-matrix(upper,nrow=15)


##Get the x and y limits (TP.FLOW, FLOW,TP)
ylim.TP<-c(min(c(mean.gamma[,3],lower[,3],upper[,3])),max(c(mean.gamma[,3],lower[,3],upper[,3])))
ylim.CROSS<-c(min(c(mean.gamma[,1],lower[,1],upper[,1])),max(c(mean.gamma[,1],lower[,1],upper[,1])))
ylim.FLOW<-c(min(c(mean.gamma[,2],lower[,2],upper[,2])),max(c(mean.gamma[,2],lower[,2],upper[,2])))
####How to make the plots in the right format
#function for th model -- #nugget +(1-exp(lag/range))
cross.sill<-model$gstat$model[[1]][2,2]
tp.sill<-model$gstat$model[[2]][2,2]
flow.sill<-model$gstat$model[[3]][2,2]
range=model$gstat$model[[1]][2,3]
x=seq(0,200,0.2)

##Make the graphics window 2 by 2
par(mfrow=c(2,2))
#First plot - the tp variogram
plot(dist,mean.gamma[,3],ylim=ylim.TP,ylab="",xlab="") #Plot the mean gamma for this variogram
lines(x,tp.sill+{1-exp(-x/range)},col="red")#Plot the actual lmcr model.
lines(dist,lower[,3],lty=2)
lines(dist,upper[,3],lty=2)
#Start the next plot at row 2 col 1
par(mfg=c(2,1))
#Plot the cross correlation variogram
plot(dist,mean.gamma[,1],ylim=ylim.TP,ylab="",xlab="") #Plot the mean gamma for this variogram
lines(x,cross.sill+{1-exp(-x/range)},col="red")#Plot the actual lmcr model.
lines(dist,lower[,1],lty=2)
lines(dist,upper[,1],lty=2)
#Plot the variogram for flow
plot(dist,mean.gamma[,2],ylim=ylim.TP,ylab="",xlab="") #Plot the mean gamma for this variogram
lines(x,flow.sill+{1-exp(-x/range)},col="red")#Plot the actual lmcr model.
lines(dist,lower[,2],lty=2)
lines(dist,upper[,2],lty=2)



###Tricking gstat into making the confidence intervals show up
###This first part will trick gstat into plotting the CI's using the directional capablilties.
##dummy.1 is the mean of the sampled variograms
dummy.1<-small.gstat[[1]]
dummy.1$dir.hor<-"mean"
dummy.1$gamma=mean.gamma
##dummy.2 is the lower CI
dummy.2<-small.gstat[[1]]
dummy.2$dir.hor<-"lower CI"
dummy.2$gamma=lower
dummy.3<-small.gstat[[1]]
dummy.3$dir.hor<-"upper CI"
dummy.3$gamma=upper
dummy<-rbind(dummy.1,dummy.2,dummy.3)
dummy$id=rep(c(rep("Total Phosphorus",15),rep("Discharge",15),rep("Total Phosphorus",15)),3)


plot(small.gstat[[1]],ids=c("one","two","there"))

 plot(rbind(small.gstat[[1]],small.gstat[[2]]),group.id=F,lty=2,pch="",col=c("blue","red"),)


 plot(rbind(small.gstat[[1]],small.gstat[[2]]),group.id=F,lty=2,pch="",col=c("blue","red"),auto.key = TRUE) 
