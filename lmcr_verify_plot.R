


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
test1<-do.call(cbind,small.gstat[-104])
##Get the gamma values from the list
gamma<-test1[,seq(3,length(test1),by=6)]

library(gstat)
##Get the mean of the gamma values
mean.gamma<-apply(gamma,1,mean)
########################################
###Plotting time. ######################
########################################
###Below is the make up for my on function to plot lmcr's. 



##And now to make a nice plotting function.


##Now put the CI and the mean values into a matrix for each variogram type order=(TP.FLOW, FLOW,TP) - this is if the gstat tricking dosent work.
mean.gamma<-matrix(mean.gamma,nrow=15)
##Get the x and y limits (TP.FLOW, FLOW,TP)
ylim.TP<-range(mean.gamma[,3])
ylim.CROSS<-range(mean.gamma[,1])
ylim.FLOW<-range(mean.gamma[,2])
####How to make the plots in the right format
#function for th model -- #nugget +(1-exp(lag/range))
cross.sill<-model$gstat$model[[1]][2,2]
tp.sill<-model$gstat$model[[2]][2,2]
flow.sill<-model$gstat$model[[3]][2,2]
cross.nugget<-model$gstat$model[[1]][1,2]
tp.nugget<-model$gstat$model[[2]][1,2]
flow.nugget<-model$gstat$model[[3]][1,2]
range=model$gstat$model[[1]][2,3]
x=seq(0,200,0.2)


##Make the graphics window 2 by 2
par(mfrow=c(2,2))
#First plot - the tp variogram
plot(dist,mean.gamma[,3],ylim=ylim.TP,ylab="semivariance",xlab="time (days)") #Plot the mean gamma for this variogram
lines(x,{tp.sill-tp.nugget}*{1-exp(-x/range)}+tp.nugget,col="red") #Plot the actual model
text(10,ylim.TP[2]*.95,"(a)")#Add the figure index
#Start the next plot at row 2 col 1
par(mfg=c(2,1))
#Plot the cross correlation variogram
plot(dist,mean.gamma[,1],ylim=ylim.CROSS,ylab="semivariance",xlab="time (days)") #Plot the mean gamma for this variogram
lines(x,{cross.sill-cross.nugget}*{1-exp(-x/range)}+cross.nugget,col="red")#Plot the actual model
text(10,ylim.CROSS[2]*.95,"(b)")#Add the figure index
#Plot the variogram for flow
plot(dist,mean.gamma[,2],ylim=ylim.FLOW,ylab="semivariance",xlab="time (days)") #Plot the mean gamma for this variogram
lines(x,{flow.sill-flow.nugget}*{1-exp(-x/range)}+flow.nugget,col="red")#Plot the actual model
text(10,ylim.FLOW[2]*.95,"(c)")#Add the figure index


