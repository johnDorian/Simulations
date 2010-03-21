####Getting the variables right.
index=0
nolags=nlags
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





### Standard deviations of each variable.
	sd[1]=sqrt(covar[1,1])
	sd[2]=sqrt(covar[2,2])

### Adjust 'AMAX' if necessary.
	amax=maxdist
	if(modtyp==1||modtyp==3||modtyp==4)amax=amax/3
	amax=amax*1.5



###	Enter initial guesses of the variogram parameters.

###	Distance parameters.
###	Nugget of the first variable's auto-variogram.
	c[1,1,1]=covar[1,1]*0.25

###	Nugget of the second variable's auto-variogram.
	c[1,2,2]=covar[2,2]*0.25

###	Nugget of the cross-variogram.
	c[1,1,2]=covar[1,2]*0.05

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


max=maximum range from fit.lmc
min=minimum range from fit.lmc
mid=max/min
perform two gstat objects and fit.lmc objects using the max and mid as ranges.
if(wss(max.fit)==wss(min.fit))return(max.fit)
if(wss(max.fit)>wss(min.fit))
	min=mid
else
	max=mid



####NOw to get this thing working. The plan is to keep changing the range with gstat.
g = gstat(g,id=c("TP","FLOW"),model=vgm(0.5,"Exp",48,0.2),fill.all=TRUE)
fit=fit.lmc(v,g,fit.lmc=FALSE,fit.range=TRUE)


max and min

#fit= gives mean and max range values to cycle through



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
			
			if(nstr>=1){
				c2=c[3,ir,ic]
				a2=a[2]
			}else{
				c2=0.0
				a2=0.0
			}
			for(j in 1:nlg){
				h=dataf[ir,ic,1,j]
				gam=dataf[ir,ic,2,j]
				rnp=dataf[ir,ic,3,j]
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
