################################################################################
### Title: lmcr_convert
### Aim: The aim of this script is to perform lmcr of two variables through time. 
### Author: Jason Lessels
### Date Created: 12/03/2010
### Last Modified: 17/03/2010
### References: The code is completly based on the work by Thomas Bishop, and is converted from the code of Lark, based on Lark's 2003 paper.
### Notes: This script has been translated from fortran code, with the prupose to make a generic function to fit a cross-co variogram with variograms supplied, by the user.
################################################################################
### WARNINGS: NO TESTING, NO IDEA IF IT WORKS...
### TODO: finish off translation, not much more to do. test to see what is wrong with it. Also have to get nvar to automatically determine how many variables exist - but might leave to down the track once everything is up and running. Get the gammah function in here. Another dream would be to get this function to accept the gstat formated cross covariogram. Create a list to wright all the results to here is an example of how to do it.
###list(WSS=1,fittingParameters=list(init.temp=1,cool.par=3),wei=3)
################################################################################
###Set the working directory (same as the git repo)
setwd("~/Documents/code/Simulations")


############################
### VARIABLE DECLARATION ###
############################

###nvar - no. variables,
###nstr - number of structures (excluding the nugget). 
###amax - max. value of distance parameter.
###icvp - set to 1 to constrain covariances to be positive, else 0.
###a - holds distance parameters.
###c - holds variances c(0:2 [structures],nvar,nvar).
###am - holds maximum change of distance parameters.
###dm - holds maximum change for each parameter (isomorphic to c).
###DATAF is (nvar,nvar, 3 [h,gamma,npairs,st. dev of gamma at h], nlag [lines]).

###Declare all the variables
lmcr<-function(semvar,nolags,nvar,wgt,icvp,cparf,modtyp,covar,maxdist,guessa,lock,c_out,a_out,vAIC)


### Final values to delcare - checked off values...
	sd=matrix(ncol=2)
	dataf=array(0,c(nvar,nvar,3,500))
###Other variables - still unsure about...
	am=matrix(ncol=2)
	c=array(NA,c(3,nvar,nvar))
	dm=array(NA,c(3,nvar,nvar))
###declaration finished...

### first loop 
	index=0
	for(ir in 1:nvar){
		for(ic in ir:nvar){
			index=index+1
			for(ix in 1:nolags){
				dataf[ir,ic,1,ix]=semvar[ix,1]
		
				if(index==1)dataf[ir,ic,2,ix]=semvar[ix,2] 
					else
				if(index==2)dataf[ir,ic,2,ix]=semvar[ix,3]
					else
				if(index==3)dataf[ir,ic,2,ix]=semvar[ix,4]
					else
				if(index==4)dataf[ir,ic,2,ix]=semvar[ix,5]
			
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
	c(1,1,2)=covar[1,2]*0.05

	if(c[1,1,1]==0||c[1,2,2]==0)c[1,1,2]=0

	c[1,2,1]=c[1,1,2]


	if(modtyp<4){
		if(modtyp==3)nstr=2 else nstr=1
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
		nstr=2
### C1 of the first variable's auto-variogram.
		c[2,1,1]=covar[1,1]*0.375

###		C1 of the second variable's auto-variogram.
		c[2,2,2]=covar[2,2]*0.375

###		C1 of the cross-variogram
		c[2,1,2]=covar[1,2]*0.1

		if(c[2,1,1]==0||c[2,2,2]==0)c[2,1,2]=0
		
		c[2,2,1]=c[1,1,2]


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
			message("")
		}	
	}
###	Enter the annealing parameters.
###	These constrain how much the variogram parameters can change with each iteration.
	for(i in 1:nstr){
		if(modtyp==1||modtyp==4){
			amin=dataf[1,1,1,2]/3
			if(a[i]<33)am[i]=1 else am[i]=2
		}
		if(modtyp==2||modtyp==5){
			amin=dataf[1,1,1,2]
			if(a[i]<100)am[i]=1 else am[i]=2
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



###	Parameter 'cparf' is the initial temperature in the input file.  If no 
###	previous attempt at a solution has been made then 'cpar' is set to 'cparf'.
	if(ics==0)cpar=cparf

	alp=0.975
	nmarkov=60
	istop=50
	iwopt=WGT
######################################################
	f<-fcn(modtyp,nlags,dataf,sd,nvar,a,c,iwopt)
	message('Weighted sum-of-squares of guessed parameters is:',f)
	message("")
###############TODO:TODO:TODO:TODO####################
###	Initiate the fitting criterion.
	call fcn(modtyp,nlags,dataf,sd,nstr,nvar,a,c,iwopt,f)
	print*, 'Weighted sum-of-squares of guessed parameters is:',f
	PRINT*,' '
	WRITE(13,*),'WSS of guessed parameters is:',f
	write(13,*) ' '
	write(13,*) 'Fitting parameters'
	write(13,*) ' '
	write(13,*) 'Initial temperature:',cpar
	write(13,*) 'Cooling parameter:',alp
	write(13,*) 'Number of trials in a Markov chain:',nmarkov
	write(13,*) 'Number of M. chains returning no change to stop:'	
     :,istop
55	continue
	write(13,*) ' '
	write(13,*) 'Weighting option:',iwopt
	write(13,*) ' '
	write(13,*) ' '
###############TODO:TODO:TODO:TODO####################

###	**************************************************************************************************

###		Start simulated annealing

###	**************************************************************************************************

###	'fic' and 'nunch' are used to count the number of chains for which 'f'
###	has not been changed.

	fic=-10.0
	nunch=0

###	Iterate cooling step.

	for(ics in 1:2000){
		if(nunch>=istop){
			break
			temp=1
		}

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
					r=rnorm(1);while(r>1||r<0){r=rnorm(1)};r
					rc=(r-0.5)*2.0*am(istr)

###					Reject inappropriate values.
					if(modtyp<4){
						if(modtyp==3&&istr==2){
							if(ao+rc)<=1||(ao+rc)>=2)next
						}else{
							if(ao+rc)<=amin||(ao+rc)>amax)next
						}
					}else{
					if(modtyp==4||modtyp==5){
						if(istr==1){
							if(ao+rc)<=amin||(ao+rc)>=a[2])next
						}else{
							if(ao+rc)<=a[1]||(ao+rc)>amax)next
						}
					}else
						break
					}
								

					a[istr]=ao+rc
					fo=f

					f<-fcn(modtyp,nlags,dataf,sd,nstr,nvar,a,c,iwopt)
					
					if(f>fo){
					
						accrej<-metrop(f,fo,cpar,accrej)
					
						if(accrej<0){
							rej=rej+1
							a[istr]=ao
							f=fo
							break
						}
					}if(accrej<0){break}
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
								if(icvp==1&&(cold+rc<0)next
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
								c[istr+1,ic,ir]=c[istr,ir,ic]
								next
							}
						}

						if(lock==1)fo=f

						f<-function(modtyp,nlags,dataf,sd,nvar,a,c,iwopt)
						for(i in 1:2){						
							if(f<=fo)break


							accrej<-metrop(f,fo,cpar)
						
							if(accrej<0){
								rej=rej+1
								c[istr+1,ir,ic]=cold
								c[istr+1,ic,ir]=c[istr+1,ir,ic]
								f=fo
								temp=2
								break
							}
						}
						if(temp==2)break					
						acc=acc+1

					}
				}
			}			

			pacc=acc/(acc+rej)

			if(ics==1)paccin=pacc
		}
		for(i in 1:2){
			if(temp==1)break
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
				iconto<-as.numeric(readline('To accept press 0, to exit press 1:'))
				if(iconto==0){
					message("")
					message("Minimising weighted sum-of-squares")
					message("(please be patient)...")
					message("")
				}else{
					if(iconto==1)stop(options("show.error.messages"=FALSE))
				}
			}

		
			cpar=cpar*alp

###		Parameter 'ics' is the Markov chain number, 'f' is the criterion minimised
###		and 'pacc' is how much it has changed since the last time.
			write(14,*)ics,f,pacc

		}

	
###	Write out results.
2001	print*,'Solution at:'
	write(13,*) 'Solution after ',ics-1,' Markov chains at:'
	
	Do 100 istr=0,nstr
		print*, 'Structure ',istr,' (0 denotes the nugget)'
		write(13,*) 'Structure ',istr,' (0 denotes the nugget)'

		if(istr.gt.0)then
			print*, 'Distance parameter:',a(istr)
			write(13,*) 'Distance parameter:',a(istr)
			A_OUT(ISTR)=DBLE(A(ISTR))
		endif

		do 110 ir=1,nvar
			do 120 ic=ir,nvar
				print*,ir,ic,c(istr,ir,ic)
				write(13,*)ir,ic,c(istr,ir,ic)
				C_OUT(istr,ir,ic)=DBLE(c(istr,ir,ic))
120			continue
110		continue
	
		IF(ISTR==2)EXIT

		print*, ' '
		print*, ' '
		write(13,*) ' '
		write(13,*) ' '

100	continue

	print*, ' '
	write(13,*) ' '

	print*,'Weighted sum-of-squares of final solution is:',f
	print*,'************************************************'
	print*,' '

	write(13,*)'Final WSS:',f


C	Overall goodness-of-fit according to Akaike Information Criterion
C	(Webster & McBratney, 1989).
	IF(MODTYP==1.OR.MODTYP==2)THEN
		vAIC=SUM(NLAGS)*LOG(F)+14
	ENDIF
	IF(MODTYP==3)THEN
		vAIC=SUM(NLAGS)*LOG(F)+16
	ENDIF
	IF(MODTYP==4.OR.MODTYP==5)THEN
		vAIC=SUM(NLAGS)*LOG(F)+22
	ENDIF
	
	WRITE(13,*)'Variable portion of Akaike Information Criterion:',vAIC
	WRITE(13,*)' '


C	Effective ranges (if applicable).
	IF(MODTYP==1.OR.MODTYP==3.OR.MODTYP==4)THEN

		IF(MODTYP==1)THEN
			WRITE(13,*)'Effective range (m) =',3*A(1)
		ENDIF
		IF(MODTYP==3)THEN
		WRITE(13,*)'Effective range (m) =',(3-((A(2)-1)*1.26795))*A(1)
		ENDIF
		IF(MODTYP==4)THEN
		WRITE(13,*)'Effective range of first structure (m) =',3*A(1)
		WRITE(13,*)'Effective range of second structure (m) =',3*A(2)
		ENDIF

		WRITE(13,*)' '
		
	ENDIF

	WRITE(13,*)'************************************************'	

	DATAF=0
	A=0
	C=0
	ICS=0

	RETURN

	END SUBROUTINE

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
	for(ir in 1:nvar){
		for(ic in ir:nvar){
			cm[ir,ic]=c[istr,ir,ic]
			cm[ic,ir]=c[istr,ic,ir]
		}
	}


### Get the eigenvalues of the matrix cm
	eval<-eigen(cm)$values
	
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

fcn<-function(modtyp,nlags,dataf,sd,nvar,a,c,iwopt){

	for(ir in 1:nvar){
		for(ic in ir:nvar){
			if(ir==1&&ic==1)nlg=nlags[1] else
			if(ir==1&&ic==2)nlg=nlags[2] else
			if(ir==2&&ic==2)nlg=nlags[3]
			C0=c[1,ir,ic]
			C1=c[1,ir,ic]
			a1=a[1]
			
			if(nstr>=1){
				c2=c[3,ir,ic]
				a2=a[1]
			}else{
				c2=0.0
				a2=0.0
			}
			for(lag in 1:nlg){
				h=dataf[ir,ic,1,lag]
				gam=dataf[ir,ic,2,lag]
				rnp=dataf[ir,ic,3,lag]
####	 TODO: NEED TO FIX UP THE NEXT FEW LINES.
c				Call 'GAMMAH' function.
				prgam=gammah(h,modtyp,C0,C1,a1,C2,a2)

### Weighted sums-of-squares
				if(ir==ic)wt=sd[ir] wt=(sd[ir])^(-4.0) else wt=((sd[ir])^(-4.0))+((sd[ic])^(-4.0))
						

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
metrop<-function(f,fo,c){
	pracc=exp((fo-f)/c)
	ran=rnorm(1);while(ran>1||ran<0){ran=rnorm(1)};ran
	if(ran<=pracc)accrej=1.0 else accrej=-1.0
	return(accrej)
}

