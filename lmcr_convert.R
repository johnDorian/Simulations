################################################################################
### Title: lmcr_convert
### Aim: The aim of this script is to perform lmcr of two variables through time. 
### Author: Jason Lessels
### Date Created: 12/03/2010
### Last Modified: 27/03/2010
### References: The code is completly based on the code supplied by Thomas Bishop which is a modification of code of Lark, based on Lark's 2003 paper.
### Notes: This script has been translated from fortran code, with the prupose to make a generic function to fit a cross-co variogram with variograms supplied, by the user.
################################################################################
### WARNINGS: LIMITED TESTING - RESULTS LOOK GOOD.
### TODO: Make variable checkers. There is very little error checking of variables.
################################################################################


############################
### VARIABLE DECLARATION ###
############################

###nvar - no. variables,
###nstr - number of structures (excluding the nugget). 
###amax - max. value of distance parameter.
###icvp - set to 1 to constrain covariances to be positive, else 0.
###a - holds distance parameters. - more than 1 value
###c - holds variances c(0:2 [structures],nvar,nvar).
###am - holds maximum change of distance parameters.
###dm - holds maximum change for each parameter (isomorphic to c).
###DATAF is (nvar,nvar, 3 [h,gamma,npairs,st. dev of gamma at h], nlag [lines]).


lmcr<-function(g,v,wgt,icvp,cpar,modtyp,covar,maxdist,guessa,lock,istop=50,plot.wss.change=TRUE){
	


###	Test each of the inputs of the function.
#####################################################
### So is g a gstat object with two variables
#Check to see if it is a gstat object
if(!is.list(g)&&class(g)[1]!="gstat")stop("variable g is not a gstat object")
#Check to see if there is two variables in the gstat object
if(length(g)>5)stop("you do not have 2 variables you need ",length(g)-5, " more variables for this function to work")
if(length(g)<5)stop("you do not have 2 variables you need ",length(g)-5, " more variables for this function to work")
### is v a variogram with the same variables as g
#Check that there is three variables in v
if(length(levels(v$id))!=3)stop("You do not have enough variables in v.")
#Check that the variables are in the right order as g
if(levels(v$id)[2]!=names(model$g)[4]||levels(v$id)[3]!=names(model$g)[5])stop(" v is not the same as g")
### is wgt between 1 and 4
if(wgt!=1)if(wgt!=2)if(wgt!=3)if(wgt!=4)stop("wgt is incorrect. Options are from 1,2,3,4.")
### is icvp 1 or 2
if(icvp!=1)if(icvp!=0)stop("icvp is incorrect. It must be 1 or 0.")
### is the cooling value positive and a integer
if(cpar<0)stop("your cooling parameter is incorrect it must be >0")
### modtyp between 1 and 8
if(modtyp!=1)if(modtyp!=2)if(modtyp!=3)if(modtyp!=4)if(modtyp!=5)if(modtyp!=6)if(modtyp!=7)if(modtyp!=8)stop("your modtyp value is incorrect it must be between 1 and 8")
### covar is the right structure (i.e. is 1,2 the same as 2,1 and a 2,2 matrix)
#Check the size of covar
if(nrow(covar)!=2)||ncol(covar)!=2)stop("covar is incorrect it must be a 2*2 matrix")
#Check the cross variogram parameters are the same
if(covar[1,2]!=covar[2,1])stop("your Cross variogarm parameter is incorrect")
### is maxdist positive
if(maxdist<0)stop("maxdist must be greater than 0")
### is guessa two variables and if modtyp is 5 or less do they both equal each other
if(length(guessa)<0)stop("no guessa provided.")
if(length(guessa)>2)stop("guessa is too long, it must be either 1 or 2 in length.")
if(length(guessa)<2&&modtyp<=5){
	guessa=c(guessa,guessa)
}
if(length(guessa)==2&&modtyp<=5){
	guessa=c(guessa[1],guessa[2])
}

### is lock 1 or 0
if(lock!=0)stop("lock must be either 1 or 0")
if(lock!=1)stop("lock must be either 1 or 0")
### is istop a integer greater than 0
#Check that it is an integer
if(istop!=floor(istop)||istop<0)stop("istop must be an integer greater than 0")
### is plot.. true or false




### g is a gstat object with both auto variables.
### v is a gstat package cross variogram.
### obtain samvar from lmcr.R script.
### nolags = 20
### wgt = ? - 1 2 3 or 4 based on the AIC function used.
### icvp = 1 - Should it return a positive definitive covar.
### cparf = - This is now cpar as the old function was based on numerous attemps from text files
### modtyp = 1:10 iso.exp=4 number based on the models suggested in gammah function.
### covar = 2*2 nuggets with minimal adjustments. - In the lmcr.R script
### maxdist = 200
### guessa = 48 48 - guess of the optimal distance for both autos. 
### lock = 0 - to make the distance the same for everything

### TODO declare a_out and c_out
###Create a list to write the results to and return at end of function. - have to declare them as zero as list wont work otherwise.
	

###	Declare the variables

### set nolags based on inputed semvariance data.
	semvar=semVar(v)
	nolags=matrix(c(length(semvar[,1]),length(semvar[,1]),length(semvar[,1])))
	nvar=2
	sd=matrix(ncol=2)
	dataf=array(0,c(nvar,nvar,3,nolags[1]))
	am=matrix(ncol=2)
	c=array(NA,c(3,nvar,nvar))
	dm=array(NA,c(3,nvar,nvar))
	temp=NA
	a_out=matrix(ncol=2)
	c_out=array(NA,c(3,nvar,nvar))


	
###declaration finished...

### first loop 
	index=0
	for(ir in 1:nvar){
		for(ic in ir:nvar){
			index=index+1
			for(ix in 1:length(nolags)){
				dataf[ir,ic,1,ix]=semvar[ix,1]
		
				if(index==1)dataf[ir,ic,2,ix]=semvar[ix,2] 
					else
				if(index==2)dataf[ir,ic,2,ix]=semvar[ix,3]
					else
				if(index==3)dataf[ir,ic,2,ix]=semvar[ix,4]
				dataf[ir,ic,3,ix]=semvar[ix,5]
			
			}
		}
	}

	nlags=nolags
	

### Standard deviations of each variable.
	sd[1]=sqrt(covar[1,1])
	sd[2]=sqrt(covar[2,2])

### Adjust 'AMAX' if necessary.
	amax=maxdist
	if(modtyp==1||modtyp==3||modtyp==4)amax=amax/3
	amax=amax*1.5



###	Enter initial guesses of the variogram parameters.

###	Distance parameters.
	a=guessa

###	Nugget of the first variable's auto-variogram.
	c[1,1,1]=covar[1,1]*0.25

###	Nugget of the second variable's auto-variogram.
	c[1,2,2]=covar[2,2]*0.25

###	Nugget of the cross-variogram.
	c[1,1,2]=covar[1,2]*0.05

	if(c[1,1,1]==0||c[1,2,2]==0)c[1,1,2]=0

	c[1,2,1]=c[1,1,2]


	if(modtyp<4){
		nstr=1
### C1 of the first variable's auto-variogram.		
		c[2,1,1]=covar[1,1]*0.75
### C1 of the second variable's auto-variogram
		c[2,2,2]=covar[2,2]*0.75
### C1 of the cross-variogram
		c[2,1,2]=covar[1,2]*0.25
		
		if(c[2,1,1]==0||c[2,1,1]==0)c[2,1,2]=0
		c[2,2,1]=c[2,1,2]
	}	



	if(modtyp>=4){
		if(modtyp>=6)nstr=2 else nstr=1
### C1 of the first variable's auto-variogram.
		c[2,1,1]=covar[1,1]*0.375

###		C1 of the second variable's auto-variogram.
		c[2,2,2]=covar[2,2]*0.375

###		C1 of the cross-variogram
		c[2,1,2]=covar[1,2]*0.1

		if(c[2,1,1]==0||c[2,2,2]==0)c[2,1,2]=0
		
		c[2,2,1]=c[2,1,2]


###		C2 of the first variable's auto-variogram.
		c[3,1,1]=covar[1,1]*0.375

###		C2 of the second variable's auto-variogram.
		c[3,2,2]=covar[2,2]*0.375

###		C2 of the cross-variogram.
		c[3,1,2]=covar[1,2]*0.25

		if(c[3,1,1]==0||c[3,2,2]==0)c[3,1,2]=0
		c[3,2,1]=c[3,1,2]
	}

###	Check guesses for positive-definiteness. 
	for(istr in 0:nstr){
		icpo<-check(c,nvar,istr)
		if(icpo<=0){
			message("")			
			message("Values for semivariogram are not positive-definite.")
			message("Try remaking semivariogram with different specifications.")
			
		}	
	}


###	Enter the annealing parameters.
###	These constrain how much the variogram parameters can change with each iteration.
	for(ix in 1:nstr){
		if(modtyp==1||modtyp==4){
			amin=dataf[1,1,1,2]/3
			if(a[ix]<33)am[ix]=1 else am[ix]=2
		}else
		if(modtyp==2||modtyp==5){
			amin=dataf[1,1,1,2]
			if(a[ix]<100)am[ix]=1 else am[ix]=2
		}
	}


	
	if(modtyp==3){
		amin=dataf[1,1,1,2]/(3-((a[2]-1)*1.26795))
		if(a[1]<(100/(3-((a[2]-1)*1.26795))))am[1]=1 else am[1]=2
		am[2]=0.05	
	}


	for(is in 0:nstr){
		for(ir in 1:nvar){
			for(ic in ir:nvar){
				if(ir==ic)
					dm[is+1,ir,ic]=covar[ir,ic]*0.25 
				else
					dm[is+1,ir,ic]=0.25*(sqrt(covar[ir,ir])*sqrt(covar[ic,ic]))
			}
		}
	}





	alp=0.975
	nmarkov=60
	iwopt=wgt

	f<-fcn(modtyp,nlags,dataf,sd,nstr,nvar,a,c,iwopt)
	message("")
	message('Weighted sum-of-squares of guessed parameters is:',f)
	message("")
	results<-list(WSSGuessedParameter=f,initialTemperature=cpar,coolingParameter=alp,
	numberTrialMarkovChain=nmarkov,numberMarkovChainReturningNoChange=istop,weightingOption=iwopt)

###	**************************************************************************************************

###		Start simulated annealing

###	**************************************************************************************************

###	'fic' and 'nunch' are used to count the number of chains for which 'f'
###	has not been changed.

	fic=-10.0
	nunch=0
### Variables for recording changing WSS.
	d.f<-data.frame(val=NA)
	d.ics<-data.frame(val=NA)
	d.pacc<-data.frame(val=NA)


###	Iterate cooling step.

	for(ics in 1:2000){

		if(nunch>istop){
			break
			temp=1
		}
		
#####	NO BRACKETS ARE BROKEN ABOVE THIS LINE
		rej=0
		acc=0

###		Iterate within the step.
		for(imc in 1:nmarkov){

###			Adjust each parameter in turn.

###			Distance parameter(s) first.
			if(lock==0){
				for(istr in 1:nstr){
					ao=a[istr]

					repeat{

						r=rnorm(1);while(r>1||r<0){r=rnorm(1)}
						rc=(r-0.5)*2.0*am[istr]

###					Reject inappropriate values.
						if(modtyp<4){
							if(modtyp==3&&istr==2){
								if((ao+rc)<=1||(ao+rc)>=2)next
							}else{
								if((ao+rc)<=amin||(ao+rc)>amax)next
							}
						}else{
							if(modtyp==4||modtyp==5){
	
								if(istr==1){
									
									if((ao+rc)<=amin||(ao+rc)>=a[2])next
									
								}else{
									if((ao+rc)<=a[1]||(ao+rc)>amax){
										
										next
									}
								}
							}
						}
						break
							
					}
				
					

					a[istr]=ao+rc
					fo=f

					f<-fcn(modtyp,nlags,dataf,sd,nstr,nvar,a,c,iwopt)
					for(i in 1:1){
						if(f<=fo)break
						accrej<-metrop(f,fo,cpar)
						if(accrej<0){
							rej=rej+1
							a[istr]=ao
							f=fo
							temp=2
							break
						}
					}
					if(!is.na(temp)&&temp==2)break					
					acc=acc+1	
				}	

			}


###			Now adjust variances.
			for(istr in 0:nstr){

				if(modtyp==3&&istr==2)break			
				for(ir in 1:nvar){
					for(ic in ir:nvar){
						cold=c[istr+1,ir,ic]
						repeat{
												
							r=rnorm(1);while(r>1||r<0){r=rnorm(1)}
							rc=(r-0.5)*2.0*dm[istr+1,ir,ic]

###						Reject inappropriate values.
							if(ir==ic){
								if((cold+rc)<0)next
							}else{
								if(icvp==1&&(cold+rc)<0)next
									else
								if(icvp==-1&&(cold+rc)>0)next
							}
						


###						Zero value if less than tolerance.
							if(abs(cold+rc)>=abs(covar[ir,ic]/1000)) c[istr+1,ir,ic]=cold+rc else c[istr+1,ir,ic]=0
							c[istr+1,ic,ir]=c[istr+1,ir,ic]


###						Check for positive-definiteness.
							icpo<-check(c,nvar,istr)
							if(icpo<0){
								c[istr+1,ir,ic]=cold
								c[istr+1,ic,ir]=c[istr+1,ir,ic]
								next
							}

						break
						}

						if(lock==1)fo=f

						f<-fcn(modtyp,nlags,dataf,sd,nstr,nvar,a,c,iwopt)


						for(i in 1:1){						
							if(f<=fo)break

							accrej<-metrop(f,fo,cpar)
							
							if(accrej<0){
								rej=rej+1
								c[istr+1,ir,ic]=cold
								c[istr+1,ic,ir]=c[istr+1,ir,ic]
								f=fo
								temp=3
								break
							}
						}
						if(!is.na(temp)&&temp==3)break					
						acc=acc+1

					}
				}
			}			

			pacc=acc/(acc+rej)

			if(ics==1)paccin=pacc
		}
		for(i in 1:1){
			if(!is.na(temp)&&temp==1)break
			if(fic==f){
				nunch=nunch+1 
			}else{
				nunch=0
				fic=f
			}

			if(ics==1){
				message(' ')
				message('Proportion of changes accepted in 1st chain: ',paccin)
				message(' ')
				message('(Between 0.90-0.99 is optimal. If equal to one, reject')
				message('and decrease the initial temperature. If less than 0.9,')
				message('reject and increase the initial temperature).')
				iconto<-as.numeric(readline('To accept press 0. to exit press 1'))

				if(iconto==0){
					message("")
					message("Minimising weighted sum-of-squares")
					message("(please be patient)...")
					message("")
					message("Warning: (Closing Graphics device will quit process)")
					
				}else{
					if(iconto==1)return(cat("Try Changing cpar and try again\n"))
				}
			}
		}

		
		cpar=cpar*alp

###		Parameter 'ics' is the Markov chain number, 'f' is the criterion minimised
###		and 'pacc' is how much it has changed since the last time.
		if(plot.wss.change==TRUE){
			if(ics==1)plot(ics,f,ylim=c(0,f*1.2),xlim=c(0,2000),cex=0.5,xlab="Perturbations",ylab="WSS")
			if(ics>1){
				points(ics,f,cex=0.5)
				legend("topright",paste("WSS= ",f),bg="white",box.col="white")
			}
		}
		d.f<-rbind(d.f,f)
		d.ics<-rbind(d.ics,ics)
		d.pacc<-rbind(d.pacc,pacc)
		
		
	}
###	Save the iteration process to the outputed list
	results$other.its=as.numeric(d.ics[,1])
	results$other.its=results$other.ics[-1]
	results$other.wss=as.numeric(d.f[,1])
	results$other.wss=results$other.f[-1]
	results$other.pacc=as.numeric(d.pacc[,1])
	results$other.pacc=results$other.pacc[-1]
############################################
###	Write out results
###Important note: To get the next part to work within the loop and add the results to the list need to first save as data.frame using rbind then convert to the list using as.numeric
	results$solutionAfter=ics-1
###	need to create a data.frame here
	structure=data.frame(val=NA)
	distance=data.frame(val=NA)
	d.ir=data.frame(val=NA)
	d.ic=data.frame(val=NA)
	d.c=data.frame(val=NA)
	for(istr in 0:nstr){
		
		structure=rbind(structure[,1],istr)
		
		if(istr>0){
			message("Solution at: ")
			message("Distance parameter: ",a[istr])
			distance=rbind(distance,a[istr])
			a_out[istr]=a[istr]
		}
		for(ir in 1:nvar){
			for(ic in ir:nvar){
				if(istr>0){				
					message(ir,ic,c[istr+1,ir,ic])
				}
				d.ir=rbind(d.ir,ir)
				d.ic=rbind(d.ic,ic)
				d.c=rbind(d.c,c[istr+1,ir,ic])
			}
		}
		if(istr==2)break
		
	}
###My new section to convert the data.frames to the lists.
	results$structure=as.numeric(structure[,1])
	results$structure=results$structure[-1]
	results$distance=as.numeric(distance[,1])
	results$distance=results$distance[-1]
	results$ir=as.numeric(d.ir[,1])
	results$ir=results$ir[-1]
	results$ic=as.numeric(d.ir[,1])
	results$ic=results$ic[-1]
	results$c=as.numeric(d.c[,1])
	results$c=results$c[-1]

	message("")
	message("Weighted sum-of-squares of final solution is: ",f)
	message("")
	
	results$finalWSS=f
###And now the variogram v added to the results
results$variogram=v
### Create a gstat object to make plotting results easier
g = gstat(g,id=c(names(g$data)[1],names(g$data)[2]),model=vgm(results$c[5],"Exp",results$distance,results$c[2]))
g = gstat(g,names(g$data)[1],model=vgm(results$c[6],"Exp",results$distance,results$c[3]),dummy=TRUE)
g = gstat(g,names(g$data)[2],model=vgm(results$c[4],"Exp",results$distance,results$c[1]),dummy=TRUE)
###Save the gstat object
results$gstat=g


###	Overall goodness-of-fit according to Akaike Information Criterion
###	(Webster & McBratney, 1989).
	if(modtyp==1||modtyp==2)vAIC=sum(nlags)*log(f)+14
	else
	if(modtyp==3)vAIC=sum(nlags)*log(f)+16
	else
	if(modtyp==4||modtyp==5)vAIC=sum(nlags)*log(f)+22
	results$variableAIC=vAIC

############################################	
###	Effective ranges (if applicable).
	if(modtyp==1||modtyp==3||modtyp==4){
		if(modtyp==1)
			results$effectiveRange1=3*a[1]
		else
		if(modtyp==3)
			results$effectiveRange2=3-((a[2]-1)*1.26795)*a[1]
		else
		if(modtyp==4){
			results$effectiveRange1stStructure=3*a[1]

		}
	}
	if(ics==2000)warning("Maximum iterations reached. Be sure to inspect final changes.")
	return(results)
}

######################################################################################################
######################################################################################################
######################################################################################################
############################### More sub rotuines - now functions ####################################
######################################################################################################
######################################################################################################
######################################################################################################


######################################################################################################
################################### The CHECK FUNCTION ###############################################
######################################################################################################


###c is the matrix in question. nvar is the number of variables. istr is not needed.
check<-function(c,nvar,istr){
	CM=array(NA,c(2,2))
	for(ir in 1:nvar){
		for(ic in ir:nvar){
			CM[ir,ic]=c[istr+1,ir,ic]
			CM[ic,ir]=c[istr+1,ic,ir]
		}
	}


### Get the eigenvalues of the matrix cm
	eval<-eigen(CM)$values
	
	rmin=10000.0
	for(ix in 1:nvar){
		if(rmin>=eval[ix]) rmin=eval[ix]
	}
	if(rmin<=-0.000001)icpo=-1 else icpo=1
	return(icpo)
}

######################################################################################################
################################### The FCN FUNCTION ###############################################
######################################################################################################
fcn<-function(modtyp,nlags,dataf,sd,nstr,nvar,a,c,iwopt){
	of=0
	rnpt=0
	for(ir in 1:nvar){
		for(ic in ir:nvar){
			if(ir==1&&ic==1)nlg=nlags[1] else
			if(ir==1&&ic==2)nlg=nlags[2] else
			if(ir==2&&ic==2)nlg=nlags[3]
			co=c[1,ir,ic]
			c1=c[2,ir,ic]
			a1=a[1]
			
			if(nstr>1){
				c2=c[3,ir,ic]
				a2=a[2]
			}else{
				c2=0.0
				a2=0.0
			}
			for(lag in 1:nlg){
				h=dataf[ir,ic,1,lag]
				gam=dataf[ir,ic,2,lag]
				rnp=dataf[ir,ic,3,lag]
####	 TODO: NEED TO FIX UP THE NEXT FEW LINES.
				prgam=gammah(h,modtyp,co,c1,a1,c2,a2)

###	Weighted sums-of-squares
				if(ir==ic)wt=(sd[ir])^(-4.0) else wt=((sd[ir])^(-4.0))+((sd[ic])^(-4.0))
				
				wdevs=((gam-prgam)^2.0)*wt

				if(iwopt==1) wdevs=wdevs else
				if(iwopt==2) wdevs=wdevs*rnp else
				if(iwopt==3) wdevs=wdevs*(rnp/(prgam^2.0)) else
				if(iwopt==4) wdevs=wdevs*(rnp/(h^2.0))
				
				of=of+wdevs
				}
			}
		}


	return(of)
}

######################################################################################################
################################### The metrop FUNCTION ##############################################
######################################################################################################
metrop<-function(f,fo,cpar){
	pracc=exp((fo-f)/cpar)
###	Using the line below it returns a random number between 0 and 1 based on the N(0,1)
	ran=rnorm(1);while(ran>1||ran<0){ran=rnorm(1)};ran
	if(ran<=pracc)accrej=1.0 else accrej=-1.0
	return(accrej)
}

######################################################################################################
################################### The semVar FUNCTION ##############################################
######################################################################################################
semVar<-function(object){
bins<-as.numeric(summary(object$id))[1]
auto.2<-split(object,object$id)[[3]]$gamma
auto.1<-split(object,object$id)[[2]]$gamma
cross<-split(object,object$id)[[1]]$gamma
dist<-split(object,object$id)[[1]]$dist
pairs<-split(object,object$id)[[1]]$np
semvar=data.frame(dist,auto.1,cross,auto.2,pairs)
}
######################################################################################################
################################### The gammah FUNCTION ##############################################
######################################################################################################
gammah<-function(h,modtyp,co,c1,a1,c2,a2){
	alpha=1.75
	sill=co+c1+c2

	switch(modtyp,
		iso.lin=return(co*h),
		iso.cir=
		if(h>a1)
			return(sill)
		else{
			hovera=h/a1
			angle=atan(sqrt(1-hovera^2)/(hovera+0.000001))
			h4=0.63662*hovera*sqrt(1-hovera^2)
			return(co+c1*(0.63662*angle+h4))
		},
		
		iso.sph=
		if(h>a1)
			return(sill)
		else{
			hovera=h/a1
			return(co+c1*(1.5*hovera-0.5*hovera^3))
		},
		iso.exp=return(co+c1*(1-exp(-h/a1))),
		
		iso.pen=if(h>a1)
			return(sill)
		else
			return(co+c1*(1.875*(h/a1)-1.25*(h/a1)^3+0.375*(h/a1)^5)),
		
		iso.dou=if(h>a2)
			return(sill)
		else
			if(h>a1)
				return(co+c1+c2*(1.5*(h/a2)-0.5*(h/a2)^3))
			else
				return(co+c1*(1.5*(h/a1)-0.5*(h/a1)^3)+c2*(1.5*(h/a2)-0.5*(h/a2)^3)),

		dou.exp=return(co+c1(1-exp((-h/a1)))+c2*(1-exp((-h/a2)))),
		
		stable=return(co+c1*(1-exp(-h^alpha/a1^alpha))),
		
		dou.sta=return(co+c1*(1-exp(-h^alpha/a1^alpha))+c2*(1-exp(-h/a2))),
		
		bro=return(co*h^a1)
	)
}
		
		
