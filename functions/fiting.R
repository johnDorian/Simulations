fiting<-function(object,parm="psill"){
	#if(class(object)[1]=="gstat")stop("object must be of class gstat")
	require(BB)
	if(parm=="nugget")parm=1 else parm=2

		gstat.num<-rep(NA,3)
		gstat.num[3]=as.numeric(object$model[[1]][,2][parm])
		gstat.num[1]=as.numeric(object$model[[2]][,2][parm])
		gstat.num[2]=as.numeric(object$model[[3]][,2][parm])
		
	f<-function(x,ini.pars){

		a11<-x[1]
		a22<-x[3]
		a12<-x[2]
		b<-rep(NA,3)
		b[1]<-(a11*a11+a12*a12)-ini.pars[1]
		b[2]<-(a22*a22+a12*a12)-ini.pars[2]
		b[3]<-(a11*a12+a22*a12)-ini.pars[3]
		return(b)
	}

	repeat{ 
		x0<-rnorm(3)
		temp<-BBsolve(par=x0,fn=f,ini.pars=gstat.num,quiet=T)
		if(temp$convergence==0){
			res<-matrix(c(temp$par[1],temp$par[2],temp$par[2],temp$par[3]),ncol=2)
			if(det(res)>=0&&min(res)>=0){
				return(res)
			}
		}
	}
}

